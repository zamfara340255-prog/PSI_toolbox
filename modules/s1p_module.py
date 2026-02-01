# -*- coding: utf-8 -*-
import math
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog
from pathlib import Path
from core.utils import tb, MATPLOTLIB_OK
from core.base_module import CalcModule

if MATPLOTLIB_OK:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    import numpy as np

class S1PVSWRModule(CalcModule):
    key = "s1p_vswr"
    title = "S1P Viewer (КСВН/S11)"

    def __init__(self, app):
        super().__init__(app)
        self.file_path = None
        self.var_path = tk.StringVar(value="Выберите файл...")
        
        # --- ОБНОВЛЕННЫЙ СПИСОК ПО УМОЛЧАНИЮ ---
        self.var_quantity = tk.StringVar(value="КСВН (VSWR)")
        
        self.var_x_unit = tk.StringVar(value="МГц")
        self.var_show_grid = tk.BooleanVar(value=True)
        self.var_threshold = tk.StringVar(value="2.0")
        self.var_cursor_mode = tk.BooleanVar(value=False)
        
        self.data_f_hz = []
        self.data_s11 = []
        self.fig, self.ax, self.canvas = None, None, None
        self._markers = [] 

    def toolbar_actions(self):
        return [
            ("Загрузить .s1p", "secondary", self.load_s1p),
            ("Автофокус (Bandwidth)", "success", self.auto_focus),
            ("Сброс зума", "info", self.reset_zoom),
            ("Открыть отдельно", "light", self.open_detached),
        ]

    def build_ui(self, parent):
        self.frame = tb.Frame(parent)
        pad = 10
        paned = tk.PanedWindow(self.frame, orient='horizontal')
        paned.pack(fill='both', expand=True)
        
        left = tb.Frame(paned, padding=pad)
        right = tb.Frame(paned, padding=pad)
        paned.add(left, width=280) 
        paned.add(right)

        # НАСТРОЙКИ
        box_file = tb.Labelframe(left, text="Данные", padding=10)
        box_file.pack(fill="x", pady=(0, 10))
        tb.Button(box_file, text="📂 Открыть .s1p", command=self.load_s1p, bootstyle="primary-outline").pack(fill="x", pady=5)
        tb.Label(box_file, textvariable=self.var_path, wraplength=220, font=("Arial", 8)).pack(anchor="w")

        box_set = tb.Labelframe(left, text="Параметры графика", padding=10)
        box_set.pack(fill="x")
        
        tb.Label(box_set, text="Тип данных:").pack(anchor="w")
        
        # --- ПОЛНЫЙ ИНЖЕНЕРНЫЙ НАБОР ---
        vals = [
            "КСВН (VSWR)", 
            "|S11| (lin)", 
            "S11 (dB)", 
            "Return Loss (dB)", 
            "Фаза (deg)",
            "Z Active (R, Ом)",   # Активное сопротивление
            "Z Reactive (X, Ом)"  # Реактивное сопротивление
        ]
        
        cb = tb.Combobox(box_set, textvariable=self.var_quantity, values=vals, state="readonly")
        cb.pack(fill="x", pady=5)
        # При смене типа данных перестраиваем график
        cb.bind("<<ComboboxSelected>>", lambda e: self.build_plot(auto_scale=True))
        
        tb.Label(box_set, text="Единицы частоты:").pack(anchor="w")
        tb.Combobox(box_set, textvariable=self.var_x_unit, values=["МГц", "ГГц"], state="readonly").pack(fill="x", pady=5)
        
        tb.Label(box_set, text="Порог (Линия/Автофокус):").pack(anchor="w")
        tb.Entry(box_set, textvariable=self.var_threshold).pack(fill="x", pady=5)

        tb.Checkbutton(box_set, text="Сетка", variable=self.var_show_grid, command=self.refresh_plot).pack(anchor="w", pady=5)
        
        # Инструменты
        box_tools = tb.Labelframe(left, text="Инструменты", padding=10)
        box_tools.pack(fill="x", pady=10)
        
        tb.Checkbutton(box_tools, text="Ставить курсор (ЛКМ)", variable=self.var_cursor_mode).pack(anchor="w")
        tb.Label(box_tools, text="ЛКМ: Вертик. линия\nПКМ: Удалить последнюю", foreground="gray", font=("Arial", 8)).pack(anchor="w")
        tb.Button(box_tools, text="Очистить все линии", command=self.clear_markers, bootstyle="danger-link").pack(anchor="w")

        # ГРАФИК
        if MATPLOTLIB_OK:
            self.fig = Figure(dpi=100)
            self.ax = self.fig.add_subplot(111)
            self.fig.subplots_adjust(bottom=0.15, left=0.15, right=0.95, top=0.92)
            
            self.canvas = FigureCanvasTkAgg(self.fig, master=right)
            self.canvas.get_tk_widget().pack(fill="both", expand=True)
            self.toolbar = NavigationToolbar2Tk(self.canvas, right)
            self.toolbar.update()
            
            self.canvas.mpl_connect("button_press_event", self.on_click)
        else:
            tb.Label(right, text="Matplotlib not found").pack()

        return self.frame

    def _parse_s1p(self, path):
        unit_mul = {"HZ":1, "KHZ":1e3, "MHZ":1e6, "GHZ":1e9}
        f_unit = "HZ"
        fh, sh = [], []
        try:
            with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    p = line.upper().strip().split()
                    if not p or p[0].startswith("!"): continue
                    if p[0].startswith("#"):
                        if "MHZ" in p: f_unit = "MHZ"
                        elif "GHZ" in p: f_unit = "GHZ"
                        elif "KHZ" in p: f_unit = "KHZ"
                        continue
                    try:
                        freq = float(p[0]) * unit_mul.get(f_unit, 1)
                        val_a = float(p[1])
                        val_b = float(p[2])
                        # Предполагаем DB Angle или Mag Angle. 
                        # Для Z нам нужен комплексный S11.
                        # Если 1-я колонка <= 0, скорее всего это dB.
                        mag = 10**(val_a/20.0) if val_a <= 0 else val_a
                        # Формируем комплексное число
                        s_complex = complex(mag * math.cos(math.radians(val_b)), mag * math.sin(math.radians(val_b)))
                        
                        sh.append(s_complex)
                        fh.append(freq)
                    except: pass
        except Exception: pass
        return fh, sh

    def load_s1p(self):
        path = filedialog.askopenfilename(filetypes=[("Touchstone", "*.s1p"), ("All", "*.*")])
        if not path: return
        self.file_path = path
        self.var_path.set(Path(path).name)
        try:
            self.data_f_hz, self.data_s11 = self._parse_s1p(path)
            self.build_plot(auto_scale=True)
        except Exception as e:
            messagebox.showerror("Err", str(e))

    def _get_ys(self, mode):
        ys = []
        Z0 = 50.0 # Характеристический импеданс
        
        for s in self.data_s11:
            mag = abs(s)
            
            if "VSWR" in mode: 
                if mag >= 0.999: v = 100.0
                else: v = (1+mag)/(1-mag)
                ys.append(v)
            elif "lin" in mode: ys.append(mag)
            elif "S11 (dB)" in mode: ys.append(20*math.log10(max(1e-9, mag)))
            elif "Return Loss" in mode: ys.append(-20*math.log10(max(1e-9, mag)))
            elif "Фаза" in mode: ys.append(math.degrees(math.atan2(s.imag, s.real)))
            
            # --- НОВЫЙ РАСЧЕТ ИМПЕДАНСА ---
            elif "Z Active" in mode or "Z Reactive" in mode:
                # Z = Z0 * (1 + S11) / (1 - S11)
                try:
                    num = 1 + s
                    den = 1 - s
                    if abs(den) < 1e-9: z = complex(1000, 1000) # Open circuit
                    else: z = Z0 * (num / den)
                    
                    if "Active" in mode: ys.append(z.real)
                    else: ys.append(z.imag)
                except: ys.append(0)
                
        return ys

    def refresh_plot(self):
        self.build_plot(auto_scale=False) 

    def build_plot(self, auto_scale=True, x_lims=None, y_lims=None, bandwidth_lines=None):
        if not self.data_f_hz or not MATPLOTLIB_OK: return
        
        self.ax.clear()
        
        div = 1e9 if self.var_x_unit.get() == "ГГц" else 1e6
        xs = [f/div for f in self.data_f_hz]
        mode = self.var_quantity.get()
        ys = self._get_ys(mode)

        # Основная линия
        self.ax.plot(xs, ys, label="Data", linewidth=1.5, color="#0056b3")

        # Лейблы
        self.ax.set_xlabel(f"Частота, {self.var_x_unit.get()}", fontsize=9, fontweight='bold')
        self.ax.set_ylabel(mode, fontsize=9, fontweight='bold')
        self.ax.grid(self.var_show_grid.get(), which='both', linestyle='-', alpha=0.7)
        self.ax.minorticks_on()
        self.ax.grid(True, which='minor', linestyle=':', alpha=0.4)

        # Линия порога (только для КСВН и RL имеет смысл рисовать красную черту)
        if "VSWR" in mode or "S11" in mode or "Return" in mode:
            try:
                thr = float(self.var_threshold.get().replace(",", "."))
                self.ax.axhline(thr, color='red', linestyle='--', linewidth=1.5, label=f"Порог {thr}")
            except: pass
        
        # Линия 50 Ом для импеданса
        if "Z Active" in mode:
            self.ax.axhline(50, color='green', linestyle='--', linewidth=1.5, label="50 Ом")
        if "Z Reactive" in mode or "Фаза" in mode:
            self.ax.axhline(0, color='green', linestyle='--', linewidth=1.5, label="0")

        # Вертикальные линии полосы
        if bandwidth_lines:
            f1, f2 = bandwidth_lines
            self.ax.axvline(f1, color='black', linestyle='--', alpha=0.6)
            self.ax.axvline(f2, color='black', linestyle='--', alpha=0.6)
            trans = self.ax.get_xaxis_transform()
            self.ax.text(f1, -0.05, f"{f1:.1f}", transform=trans, ha='center', va='top', fontsize=9, backgroundcolor='white')
            self.ax.text(f2, -0.05, f"{f2:.1f}", transform=trans, ha='center', va='top', fontsize=9, backgroundcolor='white')

        # Зум / Отступы
        if x_lims: self.ax.set_xlim(x_lims)
        elif auto_scale and xs:
            span = xs[-1] - xs[0]
            self.ax.set_xlim(xs[0] - span*0.05, xs[-1] + span*0.05)
        
        if y_lims: self.ax.set_ylim(y_lims)
        
        # Маркеры
        for x, text in self._markers:
            self.ax.axvline(x, color='green', linestyle='-.', alpha=0.8)
            # Ищем ближайший Y
            idx = min(range(len(xs)), key=lambda i: abs(xs[i]-x))
            y_val = ys[idx]
            self.ax.annotate(text, xy=(x, y_val), xytext=(5, 5), textcoords="offset points", 
                             bbox=dict(boxstyle="round", fc="w", alpha=0.8))

        self.fig.tight_layout()
        self.canvas.draw()

    def reset_zoom(self):
        self.build_plot(auto_scale=True)

    def auto_focus(self):
        """Ищет 'яму' на графике (Bandwidth) и зумит на неё."""
        if not self.data_f_hz: return
        try:
            thr = float(self.var_threshold.get().replace(",", "."))
            mode = self.var_quantity.get()
            ys = self._get_ys(mode)
            div = 1e9 if self.var_x_unit.get() == "ГГц" else 1e6
            xs = [f/div for f in self.data_f_hz]

            # Логика поиска "хорошего" диапазона
            is_vswr = "VSWR" in mode or "lin" in mode or ("S11" in mode and "dB" in mode and thr < 0)
            # Если это КСВН, ищем МЕНЬШЕ порога. Если RL, ищем БОЛЬШЕ порога.
            # Но S11 (dB) обычно < -10. 
            if "Return Loss" in mode:
                good_mask = [y >= thr for y in ys]
            elif "S11 (dB)" in mode:
                good_mask = [y <= thr for y in ys] # Порог типа -10
            else:
                good_mask = [y <= thr for y in ys] # VSWR < 2
            
            good_idxs = [i for i, val in enumerate(good_mask) if val]
            
            if not good_idxs:
                messagebox.showinfo("Результат", f"Антенна не согласована (нет точек лучше {thr})")
                return

            idx_start = good_idxs[0]
            idx_end = good_idxs[-1]
            f_start = xs[idx_start]
            f_end = xs[idx_end]
            
            bw = f_end - f_start
            margin = bw * 0.5 if bw > 0 else (xs[-1]-xs[0])*0.1
            new_xlim = (f_start - margin, f_end + margin)
            
            # Красивый масштаб по Y
            # Для КСВН: от 1.0 до порога + чуть-чуть
            if "VSWR" in mode:
                new_ylim = (1.0, max(2.5, thr + 0.5))
            else:
                new_ylim = None
            
            self.build_plot(x_lims=new_xlim, y_lims=new_ylim, bandwidth_lines=(f_start, f_end))
            self.app.set_result(f"Полоса: {f_start:.2f} - {f_end:.2f} {self.var_x_unit.get()} (BW: {bw:.2f})", "Найден резонанс", "success")

        except Exception as e:
            messagebox.showerror("Ошибка", str(e))

    def on_click(self, event):
        if not self.var_cursor_mode.get(): return
        if event.button == 1 and event.xdata: # ЛКМ
            # Форматируем текст маркера: "Freq: Val"
            y_fmt = f"{event.ydata:.2f}"
            txt = f"{event.xdata:.2f}\n{y_fmt}"
            self._markers.append((event.xdata, txt))
            self.refresh_plot()
        elif event.button == 3 and self._markers: # ПКМ
            self._markers.pop()
            self.refresh_plot()

    def clear_markers(self):
        self._markers = []
        self.refresh_plot()

    def open_detached(self):
        if not self.data_f_hz: return
        win = tk.Toplevel(self.frame)
        win.title("Детальный просмотр")
        win.geometry("1000x700")
        f = Figure(dpi=100)
        a = f.add_subplot(111)
        
        # Копируем график
        a.set_xlabel(self.ax.get_xlabel())
        a.set_ylabel(self.ax.get_ylabel())
        a.set_xlim(self.ax.get_xlim())
        a.set_ylim(self.ax.get_ylim())
        a.grid(True, which='both')
        
        for line in self.ax.get_lines():
            a.plot(line.get_xdata(), line.get_ydata(), 
                   color=line.get_color(), linestyle=line.get_linestyle(), label=line.get_label())
        
        # Копируем тексты (полоса пропускания)
        for txt in self.ax.texts:
            a.text(txt.get_position()[0], txt.get_position()[1], txt.get_text(), 
                   ha=txt.get_ha(), va=txt.get_va(), fontsize=txt.get_fontsize(), 
                   color=txt.get_color(), transform=a.get_xaxis_transform())

        a.legend()
        c = FigureCanvasTkAgg(f, master=win)
        c.get_tk_widget().pack(fill="both", expand=True)
        NavigationToolbar2Tk(c, win)
    
    def save_plot(self):
        if not self.fig: return
        path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG", "*.png")])
        if path: self.fig.savefig(path, dpi=300)
    
    def save_project(self): pass
    def load_project(self): pass
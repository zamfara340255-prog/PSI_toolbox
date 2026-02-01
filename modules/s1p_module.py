# -*- coding: utf-8 -*-
import math
import os
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog, messagebox

from core.utils import tb, MATPLOTLIB_OK
from core.base_module import CalcModule

if MATPLOTLIB_OK:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.ticker import MultipleLocator

class S1PVSWRModule(CalcModule):
    key = "s1p"
    title = "Анализ согласования (S1P)"  # Научное название

    def __init__(self, app):
        super().__init__(app)
        # --- Переменные данных ---
        self.file_path = ""
        self.raw_data_freq = [] 
        self.raw_data_s11 = []  
        
        # --- UI Переменные ---
        self.var_file = tk.StringVar(value="Файл не выбран")
        self.var_type = tk.StringVar(value="КСВН (VSWR)")
        self.var_f_unit = tk.StringVar(value="МГц")
        self.var_threshold = tk.StringVar(value="2.0")
        
        # Настройки осей
        self.var_f_start = tk.StringVar()
        self.var_f_end = tk.StringVar()
        self.var_step_x = tk.StringVar(value="100") 
        
        self.var_y_min = tk.StringVar(value="1.0")
        self.var_y_max = tk.StringVar(value="5.0")
        self.var_step_y = tk.StringVar(value="0.5")   

        self.var_markers = tk.BooleanVar(value=True)

        self.fig, self.ax, self.canvas = None, None, None
        self._markers_list = []

    def toolbar_actions(self):
        # Кнопку убрали, так как есть "Обзор" в боковой панели
        return []

    def build_ui(self, parent):
        self.frame = tb.Frame(parent)
        pad = 10
        
        paned = ttk.PanedWindow(self.frame, orient='horizontal')
        paned.pack(fill='both', expand=True)
        
        left = tb.Frame(paned, padding=pad)
        right = tb.Frame(paned, padding=pad)
        
        paned.add(left, weight=0) 
        paned.add(right, weight=1) 

        # --- ЛЕВАЯ ПАНЕЛЬ: ИСТОЧНИК ---
        lc = tb.Labelframe(left, text="Источник данных", padding=pad)
        lc.pack(fill="x", pady=(0, pad))
        
        tb.Label(lc, textvariable=self.var_file, bootstyle="secondary", wraplength=200).pack(fill="x", pady=(0, 5))
        tb.Button(lc, text="Обзор...", command=self.select_file, bootstyle="secondary-outline").pack(fill="x")

        # --- ЛЕВАЯ ПАНЕЛЬ: ПАРАМЕТРЫ ---
        pc = tb.Labelframe(left, text="Параметры", padding=pad)
        pc.pack(fill="x", pady=(0, pad))
        
        def add_combo(parent_fr, label, var, values):
            tb.Label(parent_fr, text=label, font=("Segoe UI", 9)).pack(anchor="w", pady=(5,0))
            cb = tb.Combobox(parent_fr, textvariable=var, values=values, state="readonly", font=("Segoe UI", 9))
            cb.pack(fill="x", pady=(0,5))
            cb.bind("<<ComboboxSelected>>", lambda e: self.update_plot())

        add_combo(pc, "Тип данных:", self.var_type, ["КСВН (VSWR)", "S11 (дБ)"])
        add_combo(pc, "Единицы частоты:", self.var_f_unit, ["Гц", "кГц", "МГц", "ГГц"])
        
        tb.Label(pc, text="Порог (VSWR):", font=("Segoe UI", 9)).pack(anchor="w", pady=(5,0))
        tb.Entry(pc, textvariable=self.var_threshold).pack(fill="x", pady=(0,5))

        # --- ЛЕВАЯ ПАНЕЛЬ: УПРАВЛЕНИЕ ГРАФИКОМ ---
        ac = tb.Labelframe(left, text="Управление графиком", padding=pad)
        ac.pack(fill="x", pady=(0, pad))
        
        def add_axis_row(parent_fr, txt, var):
            f = tb.Frame(parent_fr)
            f.pack(fill="x", pady=2)
            tb.Label(f, text=txt, width=12, font=("Segoe UI", 9)).pack(side="left")
            tb.Entry(f, textvariable=var, font=("Segoe UI", 9)).pack(side="right", expand=True, fill="x")

        add_axis_row(ac, "F min:", self.var_f_start)
        add_axis_row(ac, "F max:", self.var_f_end)
        add_axis_row(ac, "Шаг X:", self.var_step_x)
        tb.Separator(ac).pack(fill="x", pady=8)
        add_axis_row(ac, "Y min:", self.var_y_min)
        add_axis_row(ac, "Y max:", self.var_y_max)
        add_axis_row(ac, "Шаг Y:", self.var_step_y)

        tb.Checkbutton(ac, text="Маркеры (ЛКМ)", variable=self.var_markers).pack(anchor="w", pady=(10, 5))
        
        # Кнопки действий
        btn_frame = tb.Frame(ac)
        btn_frame.pack(fill="x", pady=(10, 0))
        
        tb.Button(btn_frame, text="Автомасштаб", command=self.autofocus, bootstyle="info").pack(fill="x", pady=2)
        tb.Button(btn_frame, text="Обновить график", command=self.update_plot, bootstyle="success").pack(fill="x", pady=2)
        tb.Button(btn_frame, text="Очистить маркеры", command=self.clear_markers, bootstyle="danger-outline").pack(fill="x", pady=2)

        # --- ПРАВАЯ ПАНЕЛЬ: ГРАФИК ---
        plot_cont = tb.Frame(right)
        plot_cont.pack(fill="both", expand=True)

        if MATPLOTLIB_OK:
            self.fig = Figure(dpi=100)
            self.ax = self.fig.add_subplot(111)
            # Убираем ручной adjust, полагаемся на tight_layout в конце
            
            self.canvas = FigureCanvasTkAgg(self.fig, master=plot_cont)
            self.canvas.get_tk_widget().pack(fill="both", expand=True)
            
            self.toolbar = NavigationToolbar2Tk(self.canvas, plot_cont)
            self.toolbar.update()
            
            self.canvas.mpl_connect("button_press_event", self.on_click)
        else:
            tb.Label(plot_cont, text="Matplotlib не установлен").pack()

        return self.frame

    # --- ЛОГИКА ---
    def select_file(self):
        path = filedialog.askopenfilename(filetypes=[("Touchstone Files", "*.s1p *.S1P"), ("All Files", "*.*")])
        if path:
            self.file_path = path
            self.var_file.set(os.path.basename(path))
            if self.parse_s1p(path):
                self.autofocus() 
            else:
                messagebox.showerror("Ошибка", "Не удалось прочитать файл .s1p")

    def parse_s1p(self, path):
        try:
            with open(path, 'r') as f:
                lines = f.readlines()
            
            freqs = []
            s11_list = []
            multiplier = 1.0
            fmt = "RI" 
            
            for line in lines:
                line = line.strip()
                if not line or line.startswith('!'): continue
                
                if line.startswith('#'):
                    parts = line.upper().split()
                    if 'HZ' in parts: multiplier = 1.0
                    if 'KHZ' in parts: multiplier = 1e3
                    if 'MHZ' in parts: multiplier = 1e6
                    if 'GHZ' in parts: multiplier = 1e9
                    
                    if 'RI' in parts: fmt = 'RI'
                    elif 'MA' in parts: fmt = 'MA'
                    elif 'DB' in parts: fmt = 'DB'
                    continue
                
                try:
                    vals = [float(x) for x in line.split()]
                    if len(vals) < 3: continue 
                    
                    f_hz = vals[0] * multiplier
                    val1, val2 = vals[1], vals[2]
                    
                    c_val = 0j
                    if fmt == 'RI':
                        c_val = complex(val1, val2)
                    elif fmt == 'MA':
                        rad = math.radians(val2)
                        c_val = complex(val1 * math.cos(rad), val1 * math.sin(rad))
                    elif fmt == 'DB':
                        mag = 10**(val1/20.0)
                        rad = math.radians(val2)
                        c_val = complex(mag * math.cos(rad), mag * math.sin(rad))
                        
                    freqs.append(f_hz)
                    s11_list.append(c_val)
                except: continue
            
            if not freqs: return False
            self.raw_data_freq = freqs
            self.raw_data_s11 = s11_list
            return True
        except Exception as e:
            print(f"S1P Parse Error: {e}")
            return False

    def get_display_factor(self):
        u = self.var_f_unit.get()
        if u == "Гц": return 1.0
        if u == "кГц": return 1e-3
        if u == "МГц": return 1e-6
        if u == "ГГц": return 1e-9
        return 1.0

    def _parse_num(self, var):
        try: return float(var.get().replace(",", "."))
        except: return None

    def calculate_crossings(self, xs, ys, threshold):
        crossings = []
        for i in range(len(ys) - 1):
            y1, y2 = ys[i], ys[i+1]
            if (y1 <= threshold < y2) or (y1 >= threshold > y2):
                x1, x2 = xs[i], xs[i+1]
                if y2 == y1: continue
                fraction = (threshold - y1) / (y2 - y1)
                x_cross = x1 + fraction * (x2 - x1)
                crossings.append(x_cross)
        return crossings

    def autofocus(self):
        if not self.raw_data_freq: return
        
        factor = self.get_display_factor()
        is_vswr = "VSWR" in self.var_type.get()
        
        x_min = self.raw_data_freq[0] * factor
        x_max = self.raw_data_freq[-1] * factor
        
        ys = []
        for s in self.raw_data_s11:
            mag = abs(s)
            if is_vswr:
                v = 100.0 if mag >= 0.999 else (1 + mag) / (1 - mag)
            else:
                v = -100.0 if mag <= 1e-9 else 20 * math.log10(mag)
            ys.append(v)
            
        if ys:
            if is_vswr:
                min_v = min(ys)
                thr = self._parse_num(self.var_threshold) or 2.0
                max_v = max(3.0, thr * 2.5)
                if min_v > max_v: max_v = min_v + 1.0
                self.var_y_min.set("1.0")
                self.var_y_max.set(f"{max_v:.2f}")
                self.var_step_y.set("0.5")
            else:
                min_v = min(ys)
                self.var_y_min.set(f"{min_v - 5:.1f}")
                self.var_y_max.set("0")
                self.var_step_y.set("5")
        
        self.var_f_start.set(f"{x_min:.4g}")
        self.var_f_end.set(f"{x_max:.4g}")
        
        span = x_max - x_min
        step_x = span / 10.0
        if step_x > 100: step_x = round(step_x / 50) * 50
        elif step_x > 10: step_x = round(step_x / 5) * 5
        else: step_x = round(step_x, 1) or 1
        self.var_step_x.set(f"{step_x}")

        self.update_plot()

    def update_plot(self):
        if not MATPLOTLIB_OK or not self.ax: return
        if not self.raw_data_freq: return
        
        try:
            self.ax.clear()
            self._markers_list = []
            
            factor = self.get_display_factor()
            xs = [f * factor for f in self.raw_data_freq]
            ys = []
            
            mode = self.var_type.get()
            is_vswr = "VSWR" in mode
            
            for s in self.raw_data_s11:
                mag = abs(s)
                if is_vswr:
                    v = 100.0 if mag >= 0.999 else (1 + mag) / (1 - mag)
                else:
                    v = -100.0 if mag <= 1e-9 else 20 * math.log10(mag)
                ys.append(v)

            # Основная линия
            self.ax.plot(xs, ys, label=mode, linewidth=2, color="#0056b3")
            
            # Порог и маркеры частот
            thr = self._parse_num(self.var_threshold)
            if thr is not None:
                self.ax.axhline(thr, color='red', linestyle='--', linewidth=1.5, label=f"Порог {thr}")
                
                crossings = self.calculate_crossings(xs, ys, thr)
                trans = self.ax.get_xaxis_transform() # Трансформация для привязки к оси X
                
                for x_cr in crossings:
                    self.ax.axvline(x_cr, color='gray', linestyle=':', linewidth=1)
                    # Подпись частоты снизу, но ВНУТРИ графика (чуть выше оси X)
                    # y=0.02 в координатах оси (где 0 - низ графика, 1 - верх)
                    self.ax.text(x_cr, 0.02, f"{x_cr:.1f}", transform=trans,
                                 rotation=90, verticalalignment='bottom', horizontalalignment='right',
                                 fontsize=8, color='#333',
                                 bbox=dict(boxstyle="square,pad=0", fc="white", alpha=0.6, ec="none"))

            # Лимиты
            x_start = self._parse_num(self.var_f_start)
            x_end = self._parse_num(self.var_f_end)
            y_min = self._parse_num(self.var_y_min)
            y_max = self._parse_num(self.var_y_max)
            
            if x_start is not None and x_end is not None: self.ax.set_xlim(x_start, x_end)
            if y_min is not None and y_max is not None: self.ax.set_ylim(y_min, y_max)
                
            # Сетка
            self.ax.grid(True, which='major', alpha=0.7)
            self.ax.grid(True, which='minor', alpha=0.3, linestyle=':')
            self.ax.minorticks_on()
            
            sx = self._parse_num(self.var_step_x)
            sy = self._parse_num(self.var_step_y)
            if sx and sx > 0: self.ax.xaxis.set_major_locator(MultipleLocator(sx))
            if sy and sy > 0: self.ax.yaxis.set_major_locator(MultipleLocator(sy))
            
            # Подписи
            u_str = self.var_f_unit.get()
            self.ax.set_xlabel(f"Частота, {u_str}", fontsize=10)
            self.ax.set_ylabel("КСВН, отн. ед." if is_vswr else "S11, дБ", fontsize=10)
            self.ax.set_title(f"График {mode}", fontsize=11)
            
            self.fig.tight_layout()
            self.canvas.draw()
            
        except Exception as e:
            print(f"Plot Error: {e}")

    def on_click(self, event):
        if not self.var_markers.get(): return
        if event.inaxes != self.ax: return
        
        if event.button == 1: 
            x, y = event.xdata, event.ydata
            pt, = self.ax.plot(x, y, 'ro', markersize=5, zorder=10)
            txt = self.ax.annotate(
                f"{x:.1f}; {y:.2f}", xy=(x, y), xytext=(10, 10), 
                textcoords='offset points', fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", fc="#ffffe0", alpha=0.9, ec="black"),
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3")
            )
            self._markers_list.append((pt, txt))
            self.canvas.draw()
            
        elif event.button == 3:
            if self._markers_list:
                pt, txt = self._markers_list.pop()
                pt.remove(); txt.remove()
                self.canvas.draw()

    def clear_markers(self):
        for pt, txt in self._markers_list:
            pt.remove()
            txt.remove()
        self._markers_list = []
        self.canvas.draw()
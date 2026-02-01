# -*- coding: utf-8 -*-
import math
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog, messagebox, simpledialog
from pathlib import Path
from core.utils import tb, MATPLOTLIB_OK
from core.base_module import CalcModule

if MATPLOTLIB_OK:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    try: import numpy as np
    except: np = None

class S1PVSWRModule(CalcModule):
    key = "s1p_vswr"
    title = "S1P Viewer (КСВН/Z/S11)"

    def __init__(self, app):
        super().__init__(app)
        self.file_path = None
        self.var_path = tk.StringVar(value="Файл не выбран")
        self.var_quantity = tk.StringVar(value="КСВН (VSWR)")
        self.var_x_unit = tk.StringVar(value="МГц")
        self.var_threshold = tk.StringVar(value="2.0")
        self.var_show_grid = tk.BooleanVar(value=True)
        self.var_cursor_mode = tk.BooleanVar(value=False)
        self.data_f_hz, self.data_s11 = [], []
        self.fig, self.ax, self.canvas = None, None, None
        self._markers = []

    def toolbar_actions(self):
        return [
            ("Загрузить .s1p", "secondary", self.load_s1p),
            ("Автофокус", "success", self.auto_focus),
            ("Сброс зума", "info", self.reset_zoom),
            ("Открыть в окне", "light", self.open_detached),
        ]

    def build_ui(self, parent):
        self.frame = tb.Frame(parent)
        pad = 10
        paned = ttk.PanedWindow(self.frame, orient='horizontal')
        paned.pack(fill='both', expand=True)
        
        left = tb.Frame(paned, padding=pad)
        right = tb.Frame(paned, padding=pad)
        paned.add(left, weight=0); paned.add(right, weight=1)

        box_file = tb.Labelframe(left, text="Файл", padding=10)
        box_file.pack(fill="x", pady=(0, 10))
        tb.Button(box_file, text="Обзор...", command=self.load_s1p).pack(fill="x", pady=5)
        tb.Label(box_file, textvariable=self.var_path, wraplength=280).pack(anchor="w")

        box_set = tb.Labelframe(left, text="Параметры", padding=10)
        box_set.pack(fill="x", pady=(0, 10))
        
        tb.Label(box_set, text="Тип данных:").pack(anchor="w")
        vals = ["КСВН (VSWR)", "|S11| (lin)", "S11 (dB)", "Return Loss (dB)", "Фаза (deg)", "Z Active (R)", "Z Reactive (X)"]
        cb = tb.Combobox(box_set, textvariable=self.var_quantity, values=vals, state="readonly")
        cb.pack(fill="x", pady=5)
        cb.bind("<<ComboboxSelected>>", lambda e: self.build_plot(auto_scale=True))
        
        tb.Label(box_set, text="Единицы частоты:").pack(anchor="w")
        # Расширенный выбор единиц
        cb_f = tb.Combobox(box_set, textvariable=self.var_x_unit, values=["Гц", "кГц", "МГц", "ГГц"], state="readonly")
        cb_f.pack(fill="x", pady=5)
        cb_f.bind("<<ComboboxSelected>>", lambda e: self.build_plot(auto_scale=True))
        
        tb.Label(box_set, text="Порог:").pack(anchor="w")
        tb.Entry(box_set, textvariable=self.var_threshold).pack(fill="x", pady=5)
        tb.Checkbutton(box_set, text="Сетка", variable=self.var_show_grid, command=self.refresh_plot).pack(anchor="w")
        
        box_tools = tb.Labelframe(left, text="Курсор", padding=10)
        box_tools.pack(fill="x")
        tb.Checkbutton(box_tools, text="Ставить ЛКМ", variable=self.var_cursor_mode).pack(anchor="w")
        tb.Button(box_tools, text="Очистить", command=self.clear_markers, bootstyle="danger-link").pack(anchor="w")

        if MATPLOTLIB_OK:
            self.fig = Figure(dpi=100)
            self.ax = self.fig.add_subplot(111)
            # Отступы
            self.fig.subplots_adjust(left=0.12, bottom=0.12, right=0.96, top=0.93)
            self.canvas = FigureCanvasTkAgg(self.fig, master=right)
            self.canvas.get_tk_widget().pack(fill="both", expand=True)
            self.toolbar = NavigationToolbar2Tk(self.canvas, right)
            self.toolbar.update()
            self.canvas.mpl_connect("button_press_event", self.on_click)
        else: tb.Label(right, text="No Matplotlib").pack()
        return self.frame

    def _parse_s1p(self, path):
        mul = {"HZ":1, "KHZ":1e3, "MHZ":1e6, "GHZ":1e9}
        unit, fmt = "HZ", "MA"
        fh, sh = [], []
        try:
            with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    p = line.upper().strip().split()
                    if not p or p[0].startswith("!"): continue
                    if p[0].startswith("#"):
                        if "HZ" in p: unit="HZ"
                        if "KHZ" in p: unit="KHZ"
                        if "MHZ" in p: unit="MHZ"
                        if "GHZ" in p: unit="GHZ"
                        if "RI" in p: fmt="RI"
                        if "DB" in p: fmt="DB"
                        continue
                    try:
                        f_val = float(p[0]) * mul.get(unit, 1)
                        a, b = float(p[1]), float(p[2])
                        s = complex(0,0)
                        if fmt=="RI": s=complex(a,b)
                        elif fmt=="DB": 
                            mag = 10**(a/20.0)
                            s=mag*complex(math.cos(math.radians(b)), math.sin(math.radians(b)))
                        else: # MA
                            s=a*complex(math.cos(math.radians(b)), math.sin(math.radians(b)))
                        fh.append(f_val); sh.append(s)
                    except: pass
        except: pass
        return fh, sh

    def load_s1p(self):
        path = filedialog.askopenfilename(filetypes=[("S1P", "*.s1p"), ("All", "*.*")])
        if not path: return
        self.file_path = path
        self.var_path.set(Path(path).name)
        self.data_f_hz, self.data_s11 = self._parse_s1p(path)
        self.build_plot(auto_scale=True)

    def _get_ys(self, mode):
        ys = []
        Z0 = 50.0
        for s in self.data_s11:
            mag = abs(s)
            if "VSWR" in mode: 
                ys.append((1+mag)/(1-mag) if mag<0.999 else 100)
            elif "lin" in mode: ys.append(mag)
            elif "dB" in mode and "S11" in mode: ys.append(20*math.log10(max(1e-9, mag)))
            elif "Return" in mode: ys.append(-20*math.log10(max(1e-9, mag)))
            elif "Фаза" in mode: ys.append(math.degrees(math.atan2(s.imag, s.real)))
            elif "Active" in mode:
                try: z=Z0*((1+s)/(1-s)); ys.append(z.real)
                except: ys.append(0)
            elif "Reactive" in mode:
                try: z=Z0*((1+s)/(1-s)); ys.append(z.imag)
                except: ys.append(0)
        return ys

    def refresh_plot(self): self.build_plot(auto_scale=False)
    def reset_zoom(self): self.build_plot(auto_scale=True)

    def build_plot(self, auto_scale=True, x_lims=None, bandwidth_lines=None):
        if not self.data_f_hz or not MATPLOTLIB_OK: return
        self.ax.clear()
        
        # Пересчет единиц
        u = self.var_x_unit.get()
        div = 1e9 if u=="ГГц" else (1e6 if u=="МГц" else (1e3 if u=="кГц" else 1))
        
        xs = [f/div for f in self.data_f_hz]
        mode = self.var_quantity.get()
        ys = self._get_ys(mode)
        
        self.ax.plot(xs, ys, label="Data", linewidth=1.5, color="#0056b3")
        self.ax.set_xlabel(f"Частота ({u})", fontsize=9, fontweight='bold')
        self.ax.set_ylabel(mode, fontsize=9, fontweight='bold')
        self.ax.grid(self.var_show_grid.get())

        # Порог
        try:
            thr = float(self.var_threshold.get().replace(",", "."))
            if "VSWR" in mode: self.ax.axhline(thr, color='red', linestyle='--')
        except: pass

        if x_lims: self.ax.set_xlim(x_lims)
        elif auto_scale and xs:
            span = max(xs)-min(xs)
            m = span*0.05 if span>0 else 1
            self.ax.set_xlim(min(xs)-m, max(xs)+m)

        # Линии полосы
        if bandwidth_lines:
            f1, f2 = bandwidth_lines
            # Делим f1, f2 на div, так как они в Гц
            self.ax.axvline(f1/div, color='black', linestyle='--')
            self.ax.axvline(f2/div, color='black', linestyle='--')
            # Подписи
            trans = self.ax.get_xaxis_transform()
            self.ax.text(f1/div, -0.05, f"{f1/div:.2f}", transform=trans, ha='center', va='top', fontsize=9, backgroundcolor='white')
            self.ax.text(f2/div, -0.05, f"{f2/div:.2f}", transform=trans, ha='center', va='top', fontsize=9, backgroundcolor='white')

        for mx, mt in self._markers:
            self.ax.axvline(mx, color='green', linestyle=':')
            if len(xs)>1:
                idx = min(range(len(xs)), key=lambda i: abs(xs[i]-mx))
                self.ax.annotate(mt, xy=(mx, ys[idx]), xytext=(5,5), textcoords="offset points", bbox=dict(boxstyle="round", fc="white"))
        
        self.canvas.draw()

    def auto_focus(self):
        # Старая добрая логика: ищем точки, берем min/max, зумим
        if not self.data_f_hz: return
        try:
            thr = float(self.var_threshold.get().replace(",", "."))
            mode = self.var_quantity.get()
            ys = self._get_ys(mode)
            
            u = self.var_x_unit.get()
            div = 1e9 if u=="ГГц" else (1e6 if u=="МГц" else (1e3 if u=="кГц" else 1))
            
            # Работаем в Гц для поиска
            xs_hz = self.data_f_hz

            # Условие: КСВН < порога
            is_good = lambda y: y <= thr if "VSWR" in mode else (y >= thr if "Return" in mode else y <= thr)
            
            good_idxs = [i for i, y in enumerate(ys) if is_good(y)]
            if not good_idxs:
                messagebox.showinfo("Инфо", f"Нет точек, удовлетворяющих условию {thr}")
                return
            
            # Берем просто минимум и максимум из хороших точек (в Гц)
            f_start_hz = xs_hz[min(good_idxs)]
            f_end_hz = xs_hz[max(good_idxs)]
            
            # Переводим в текущие единицы для графика
            f_start = f_start_hz / div
            f_end = f_end_hz / div
            
            bw = f_end - f_start
            m = bw*0.1 if bw>0 else (f_end)*0.05
            
            # Вызываем с линиями (передаем Гц, build_plot сам поделит)
            self.build_plot(auto_scale=False, x_lims=(f_start-m, f_end+m), bandwidth_lines=(f_start_hz, f_end_hz))
            
        except Exception as e: messagebox.showerror("Err", str(e))

    def on_click(self, event):
        if not self.var_cursor_mode.get(): return
        if event.button==1 and event.xdata:
            self._markers.append((event.xdata, f"{event.xdata:.2f}"))
            self.refresh_plot()
        elif event.button==3 and self._markers:
            self._markers.pop(); self.refresh_plot()
            
    def clear_markers(self): self._markers=[]; self.refresh_plot()
    def open_detached(self): pass # (как раньше)
    def save_plot(self): pass
    def save_project(self): pass
    def load_project(self): pass
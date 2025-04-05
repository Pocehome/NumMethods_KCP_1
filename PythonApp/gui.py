import tkinter as tk
from tkinter import messagebox
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from KSR1_Module import KSR_Task


class Window:
    def __init__(self, master):
        # variables
        self.table_data = None
        self.A = None
        self.B = None
        self.U0 = None
        self.step_size = None
        self.max_e_error = None
        self.e_border = None
        self.max_steps = None
        self.k = None
        self.k_with_star = None
        self.c = None
        self.m = None
        self.final_ref = None
        self.mode = 1
        self.graph_mode = 0

        # window master
        self.master = master
        self.master.title("КСР №11 вариант №3")
        self.master.geometry("1300x600")

        # Toggle button
        tk.Label(self.master, text="Режим работы:").grid(row=4, column=0, sticky='e')
        self.button_mode = tk.Button(self.master, text="С ОЛП", command=self.toggle_mode)
        self.button_mode.grid(row=4, column=1, sticky='we')

        self.create_KSR_task_widgets()

    def clear_widgets(self):
        # Remove all widgets except the dropdown
        for widget in self.master.winfo_children():
            if widget not in [self.button_mode]:
                widget.destroy()

        # Clear variables
        self.table_data = None
        self.A = None
        self.B = None
        self.U0 = None
        self.step_size = None
        self.max_e_error = None
        self.e_border = None
        self.max_steps = None
        self.k = None
        self.k_with_star = None
        self.c = None
        self.m = None
        self.final_ref = None
        self.graph_mode = 0

    def create_KSR_task_widgets(self):
        self.master.grid_columnconfigure(0, weight=200)
        self.master.grid_columnconfigure(1, weight=1)
        self.master.grid_columnconfigure(2, weight=10)
        self.master.grid_columnconfigure(3, weight=1)

        # Input fields for test task
        tk.Label(self.master, text="A:").grid(row=0, column=2, sticky='e')
        tk.Label(self.master, text="B:").grid(row=1, column=2, sticky='e')
        tk.Label(self.master, text="U0:").grid(row=2, column=2, sticky='e')
        tk.Label(self.master, text="Начальный шаг:").grid(row=3, column=2, sticky='e')
        tk.Label(self.master, text="Макс. число шагов:").grid(row=4, column=2, sticky='e')
        tk.Label(self.master, text="Точность выхода на границу:").grid(row=5, column=2, sticky='e')
        if self.mode == 1:
            tk.Label(self.master, text="Контроль лок. погрешности:").grid(row=6, column=2, sticky='e')
        tk.Label(self.master, text="k:").grid(row=0, column=0, sticky='e')
        tk.Label(self.master, text="k*:").grid(row=1, column=0, sticky='e')
        tk.Label(self.master, text="c:").grid(row=2, column=0, sticky='e')
        tk.Label(self.master, text="m:").grid(row=3, column=0, sticky='e')

        # Input fields
        self.entry_a = tk.Entry(self.master)
        self.entry_b = tk.Entry(self.master)
        self.entry_u0 = tk.Entry(self.master)
        self.entry_step_size = tk.Entry(self.master)
        self.entry_max_steps = tk.Entry(self.master)
        self.entry_e_border = tk.Entry(self.master)
        if self.mode == 1:
            self.entry_max_e_error = tk.Entry(self.master)
        self.entry_k = tk.Entry(self.master)
        self.entry_k_with_star = tk.Entry(self.master)
        self.entry_c = tk.Entry(self.master)
        self.entry_m = tk.Entry(self.master)

        # Set default values
        self.entry_a.insert(0, "0.0")
        self.entry_b.insert(0, "2.0")
        self.entry_u0.insert(0, "10.0")
        self.entry_step_size.insert(0, "0.01")
        self.entry_max_steps.insert(0, "100000")
        self.entry_e_border.insert(0, "1e-6")
        if self.mode == 1:
            self.entry_max_e_error.insert(0, "1e-12")
        self.entry_k.insert(0, "2")
        self.entry_k_with_star.insert(0, "2")
        self.entry_c.insert(0, "0.15")
        self.entry_m.insert(0, "0.01")

        # Placement of input fields
        self.entry_a.grid(row=0, column=3, sticky='ew')
        self.entry_b.grid(row=1, column=3, sticky='ew')
        self.entry_u0.grid(row=2, column=3, sticky='ew')
        self.entry_step_size.grid(row=3, column=3, sticky='ew')
        self.entry_max_steps.grid(row=4, column=3, sticky='ew')
        self.entry_e_border.grid(row=5, column=3, sticky='ew')
        if self.mode == 1:
            self.entry_max_e_error.grid(row=6, column=3, sticky='ew')
        self.entry_k.grid(row=0, column=1, sticky='ew')
        self.entry_k_with_star.grid(row=1, column=1, sticky='ew')
        self.entry_c.grid(row=2, column=1, sticky='ew')
        self.entry_m.grid(row=3, column=1, sticky='ew')

        # Calculate button
        self.calculate_button = tk.Button(self.master, text="Вычислить", command=self.calculate_values)
        self.calculate_button.grid(row=6, column=0, columnspan=2, sticky='we')

        # Create Treeview to display the table
        self.tree = ttk.Treeview(self.master, columns=["Col" + str(i) for i in range(10)], show="headings")
        self.tree.grid(row=7, rowspan=2, column=0, columnspan=4, sticky='nsew', padx=5,
                       pady=5)  # Place the table under input fields

        # Add horizontal and vertical scrolling
        vsb = ttk.Scrollbar(self.master, orient="vertical", command=self.tree.yview)
        vsb.grid(row=7, column=3, sticky='sne', rowspan=3)
        self.tree.configure(yscrollcommand=vsb.set)

        hsb = ttk.Scrollbar(self.master, orient="horizontal", command=self.tree.xview)
        hsb.grid(row=9, column=0, columnspan=4, sticky='esw')
        self.tree.configure(xscrollcommand=hsb.set)

        # Configure stretching
        self.master.grid_rowconfigure(8, weight=1)
        self.master.grid_columnconfigure(0, weight=1)

        # Create figure for the plot
        self.figure = Figure(figsize=(6, 4), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.master)

        # Create a Text widget for final report
        self.final_report_text = tk.Text(self.master, height=10, width=50)
        self.final_report_text.grid(row=0, column=4, rowspan=7, sticky='nwe', padx=5, pady=5)

        # Place the canvas (the graph) in the bottom right
        self.canvas.get_tk_widget().grid(row=8, column=4, rowspan=5, sticky='sewn', padx=5, pady=5)

        # Кнопка для построения графика
        self.plot_button = tk.Button(self.master, text="Построить график V2(x)", command=self.toggle_plot)
        self.plot_button.grid(row=7, column=4, sticky='we', padx=5, pady=5)

    def calculate_values(self):
        try:
            # Получение значений из полей ввода
            self.A = float(self.entry_a.get())
            self.B = float(self.entry_b.get())
            self.U0 = float(self.entry_u0.get())
            self.step_size = float(self.entry_step_size.get())
            self.max_steps = int(self.entry_max_steps.get())
            self.e_border = float(self.entry_e_border.get())
            if self.mode == 1:
                self.max_e_error = float(self.entry_max_e_error.get())
            else:
                self.max_e_error = 0

            self.k = float(self.entry_k.get())
            self.k_with_star = float(self.entry_k_with_star.get())
            self.c = float(self.entry_c.get())
            self.m = float(self.entry_m.get())

            # Проверки на корректность значений
            if self.A >= self.B:
                messagebox.showerror("Ошибка", "Значение A должно быть меньше B.")
                return
            if self.max_steps <= 1:
                messagebox.showerror("Ошибка", "Максимальное количество шагов должно быть больше 1.")
                return

            task = KSR_Task(self.A, self.B, self.step_size, self.max_e_error, self.e_border, self.max_steps, self.U0)
            task.set_params(self.k, self.k_with_star, self.c, self.m)

            # Вычисление
            if self.mode == 1:
                task.Solve_With_Error_Control()
            else:
                task.Solve_Without_Error_Control()
            self.table_data = task.get_table_information()
            self.update_table()

            # Отчет
            self.final_ref = task.get_final_reference()
            self.show_final_reference()

            # Обновление графика
            if self.graph_mode == 0:
                self.plot_graph_V1_x()
            elif self.graph_mode == 1:
                self.plot_graph_V2_x()
            else:
                self.plot_graph_V2_V1()

            # print(task.get_ERRORS_LIST())

        except ValueError:
            messagebox.showerror("Ошибка", "Пожалуйста, введите корректные значения.")

    # Функция для переключения между режимами
    def toggle_mode(self):
        self.clear_widgets()

        # Toggle button
        tk.Label(self.master, text="Режим работы:").grid(row=4, column=0, sticky='e')
        self.button_mode = tk.Button(self.master, command=self.toggle_mode)
        self.button_mode.grid(row=4, column=1, sticky='we')

        if self.mode == 1:
            self.button_mode.config(text='Без ОЛП')
            self.mode = 0
        else:
            self.button_mode.config(text='С ОЛП')
            self.mode = 1

        self.create_KSR_task_widgets()

    # Функция для переключения между графиками
    def toggle_plot(self):
        if self.graph_mode == 0:
            self.plot_graph_V2_x()
            self.plot_button.config(text="Построить график V2(V1)")
            self.graph_mode = 1

        elif self.graph_mode == 1:
            self.plot_graph_V2_V1()
            self.plot_button.config(text="Построить график V1(x)")
            self.graph_mode = 2

        elif self.graph_mode == 2:
            self.plot_graph_V1_x()
            self.plot_button.config(text="Построить график V2(x)")
            self.graph_mode = 0

    def update_table(self):
        # Clear Treeview before updating
        for row in self.tree.get_children():
            self.tree.delete(row)

        # Set column headers
        if self.mode == 0:
            columns = ["i", "X", "V1", "V2", "h"]
            arr_width = [50, 100, 100, 100, 100]
        else:
            columns = ["i", "X", "V1", "V2", "V1^", "V2^", "V1-V1^", "V2-V2^", "ОЛП", "h", "/2", "*2"]
            arr_width = [50, 100, 100, 100, 100, 100, 100, 100, 100, 100, 50, 50]

        if list(self.tree["columns"]) != columns:
            self.tree["columns"] = columns
            for i in range(len(columns)):
                self.tree.heading(columns[i], text=columns[i])
                self.tree.column(columns[i], width=arr_width[i], minwidth=arr_width[i])  # Set column width

        # Add data to the table with formatting
        for row in self.table_data:
            formatted_row = [f"{value:.6g}" if isinstance(value, float) else value for value in row]
            self.tree.insert("", "end", values=formatted_row)

    def plot_graph_V1_x(self):
        # Очистка графика
        self.ax.clear()

        # Extract X and V values from table data for the second plot
        X_values = [row[1] for row in self.table_data]
        V1_values = [row[2] for row in self.table_data]

        # Построение графика численного решения
        self.ax.plot(X_values, V1_values, label='V1(x)', color='red', alpha=0.7)

        # Customize the plot
        self.ax.set_xlabel('x')
        # self.ax.xaxis.set_tick_params(labelsize=8)
        # self.ax.yaxis.set_tick_params(labelsize=8)
        self.ax.legend()
        self.ax.grid()

        # Автоматическая настройка размещения элементов графика
        self.ax.figure.tight_layout()

        # Update the plot
        self.canvas.draw()

    def plot_graph_V2_x(self):
        # Очистка графика
        self.ax.clear()

        # Extract X and V values from table data for the second plot
        X_values = [row[1] for row in self.table_data]
        V2_values = [row[3] for row in self.table_data]

        # Построение графика численного решения
        self.ax.plot(X_values, V2_values, label='V2(x)', color='red', alpha=0.7)

        # Customize the plot
        self.ax.set_xlabel('x')
        # self.ax.xaxis.set_tick_params(labelsize=8)
        # self.ax.yaxis.set_tick_params(labelsize=8)
        self.ax.legend()
        self.ax.grid()

        # Автоматическая настройка размещения элементов графика
        self.ax.figure.tight_layout()

        # Update the plot
        self.canvas.draw()

    def plot_graph_V2_V1(self):
        try:
            # Extract X and V values from table data for the second plot
            V1_values = [row[2] for row in self.table_data]
            V2_values = [row[3] for row in self.table_data]

            # Очистка и построение графика
            self.ax.clear()
            self.ax.plot(V1_values, V2_values, label="V2(V1)", color='blue')

            # Customize the plot
            self.ax.set_xlabel("V1")
            # self.ax.set_ylabel("V2")
            self.ax.legend()
            self.ax.grid()

            # Update the plot
            self.canvas.draw()

        except Exception as e:
            messagebox.showerror("Ошибка", f"Не удалось построить график: {e}")

    def show_final_reference(self):
        # Clear the previous report text
        self.final_report_text.delete(1.0, tk.END)

        # Insert new information into the Text widget
        if self.mode == 0:
            info = (
                f"Кол-во итераций: {self.final_ref.ITERATIONS_COUNT}\n"
                f"Расстояние до последней точки: {self.final_ref.DISTANCE_B_LAST_POINT}\n"
            )
        else:
            info = (
                f"Кол-во итераций: {self.final_ref.ITERATIONS_COUNT}\n"
                f"Расстояние до последней точки: {self.final_ref.DISTANCE_B_LAST_POINT}\n"
                f"Кол-во удвоений шага: {self.final_ref.STEP_DOUBLING_COUNT}\n"
                f"Кол-во уменьшений шага: {self.final_ref.STEP_REDUCTION_COUNT}\n"
                f"Макс. шаг с X: {self.final_ref.MAX_STEP_WITH_X}\n"
                f"Мин. шаг с X: {self.final_ref.MIN_STEP_WITH_X}\n"
                f"Последнее значение X: {self.final_ref.LAST_X}\n"
                f"Последнее значение V: {self.final_ref.LAST_V}\n"
                f"Макс. значение оценки олп: {self.final_ref.MAX_ERROR}\n"
                f"Мин. значение оценки олп: {self.final_ref.MIN_ERROR}\n"

            )
        if self.final_ref.IS_INF:
            info += (
                "На следующем шаге значение X перейдёт через асимптоту, поэтому вычисление следующего значения V невозможно"
            )

        self.final_report_text.insert(tk.END, info)


def create_gui():
    root = tk.Tk()
    app = Window(root)
    root.mainloop()


if __name__ == "__main__":
    create_gui()
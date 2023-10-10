import time

from gui import Ui_Dialog
from data_and_processing import DataAndProcessing
from graph import Graph
from  drawer import Drawer

# КЛАСС АЛГОРИТМА ПРИЛОЖЕНИЯ
class GuiProgram(Ui_Dialog):

    def __init__(self, dialog):
        # ПОЛЯ КЛАССА
        self.data_and_processing = DataAndProcessing()

        # ДЕЙСТВИЯ ПРИ ВКЛЮЧЕНИИ
        # Создаем окно
        Ui_Dialog.__init__(self)
        self.setupUi(dialog)  # Устанавливаем пользовательский интерфейс

        # Параметры 1 графика
        self.graph_1 = Graph(
            layout=self.layout_plot_1,
            widget=self.widget_plot_1
        )

        self.pushButton_initial_state.clicked.connect(self.initial_state)
        self.pushButton_start_processing.clicked.connect(self.start_processing)

    def initial_state(self):
        start_time = time.perf_counter()
        self.data_and_processing.initial_state(
            number_elements_on_side=int(self.lineEdit_elements_on_side.text()),
            start_x=float(self.lineEdit_x0.text()),
            start_y=float(self.lineEdit_y0.text()),
            temperature=float(self.lineEdit_temperature.text()) # Кельвинов
        )
        end_time = time.perf_counter()
        print(end_time - start_time, "seconds")
        Drawer.particle_field(self.graph_1, self.data_and_processing)

    def start_processing(self):
        start_time = time.perf_counter()
        self.data_and_processing.state_calculation(
            time_step=int(self.lineEdit_time_step.text()),
            number_of_steps=int(self.lineEdit_iterations.text()),
        )

        end_time = time.perf_counter()
        print(end_time - start_time, "seconds")
        Drawer.particle_field(self.graph_1, self.data_and_processing)

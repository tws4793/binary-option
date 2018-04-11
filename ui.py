import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_qt5agg import \
    FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5 import uic
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

from bin_option import *
from option import *
from utils import *

"""
Desktop User Interface
Using PyQt5
"""

# Constants
APP_NAME = 'Option Pricing'
WINDOW_SIZE = 1500,500
OPTION_OPTIONS = 'Call','Put'
OPTION_VARS = 'S','K', 'r', 'sigma', 'T', 'q'
OPTION_ITEMS = 'Black Scholes', 'Binomial', 'Monte-Carlo', 'Implicit', 'Explicit', 'Crank-Nicolson'
OPTION_TYPES = 'Vanilla', 'Binary Cash', 'Binary Asset'
BTN_UPDATE = 'Update'
BTN_CALCULATE = 'Calculate'
BTN_RESET = 'Reset'
SIMULATIONS = 50_000

# option_values = [0,0,0,0.0,0,0]
# option_values = [50,50,0.017,0.4,(183/365),0.01]
option_values = [0,0,0.017,0,(183/365),0]

display_type = {
    'vanilla': 0.0,
    'cash': 0.0,
    'asset': 0.0,
}

class Main(QWidget):
    def __init__(self, *args, **kwargs):
        super(Main, self).__init__(*args,**kwargs)
        self.setup_ui()
    
    def setup_ui(self):
        """Performs the setup of the UI
        """
        self.resize(WINDOW_SIZE[0], WINDOW_SIZE[1])
        self.setWindowTitle(APP_NAME)

        # Overall
        layout = QHBoxLayout()
        layout.addLayout(self.layout_input())
        layout.addLayout(self.layout_output())
        self.setLayout(layout)

    def layout_input(self):
        """This method is for the setup of the output section.

        Returns:
            The layout for the input section.

        """
        # Setup the layouts
        layout = QFormLayout()

        # == Ticker ==
        self.tb_ticker = QLineEdit(self)
        self.tb_ticker.setMaxLength(4)
        layout.addRow('Ticker', self.tb_ticker)
        
        self.btn_get_ticker_data = QPushButton(BTN_UPDATE)
        layout.addWidget(self.btn_get_ticker_data)
        self.lb_update_progress = QLabel('')
        layout.addWidget(self.lb_update_progress)
        self.btn_get_ticker_data.clicked.connect(self.get_data)
        layout.addWidget(QSplitter(Qt.Vertical))
        
        # == Option Pricing Variables ==
        self.tb_v = [self.display_input_text(str(option_values[i])) for i,tb in enumerate(OPTION_VARS)]
        for i,v in enumerate(OPTION_VARS):
            layout.addRow(v, self.tb_v[i])
        self.btn_calculate = QPushButton(BTN_CALCULATE)
        # layout.addWidget(self.btn_calculate)
        self.btn_calculate.clicked.connect(self.calculate)

        # == Reset ==
        self.btn_reset = QPushButton(BTN_RESET)
        self.btn_reset = QPushButton()
        # layout.addWidget(self.btn_reset)
        self.btn_reset.clicked.connect(self.reset_data)
        layout.addRow(self.btn_reset, self.btn_calculate)

        return layout
    
    def layout_output(self):
        """This method is for the setup of the output section.

        Returns:
            The layout for the output section.

        """
        layout = QVBoxLayout()

        # == Manipulate Output ==
        layout_manipulation = QHBoxLayout()

        # Models
        self.cb_models = QComboBox(self)
        self.cb_models.addItems(OPTION_ITEMS)
        layout_manipulation.addWidget(self.cb_models)
        
        # Options
        self.mode_cp = 0
        self.btn_cp = QPushButton(OPTION_OPTIONS[self.mode_cp])
        self.btn_cp.clicked.connect(self.toggle_cp)
        layout_manipulation.addWidget(self.btn_cp)

        # Consolidate all layouts
        layout.addLayout(layout_manipulation)

        # == Graphs ==
        layout_graphs = QHBoxLayout()

        # Prepare the layouts
        self.lb_header = [QLabel(v) for i,v in enumerate(OPTION_TYPES)]
        self.figures = [plt.figure() for i,v in enumerate(OPTION_TYPES)]
        self.graphs = [FigureCanvas(self.figures[i]) for i,v in enumerate(self.figures)]
        self.lb_figure = [QLabel('${:,.2f}'.format(v[1])) for k,v in enumerate(display_type.items())]

        # Add the individual layouts
        for i,v in enumerate(OPTION_TYPES):
            self.lb_header[i].setAlignment(Qt.AlignCenter)
            self.lb_figure[i].setAlignment(Qt.AlignCenter)
            layout_graph = QVBoxLayout()
            layout_graph.addWidget(self.lb_header[i])
            layout_graph.addWidget(self.graphs[i])
            layout_graph.addWidget(self.lb_figure[i])
            layout_graphs.addLayout(layout_graph) # Add

        # Consolidate all layouts
        layout.addLayout(layout_graphs)
        
        return layout
    
    def display_input_text(self,default_text):
        """Displays the input box
        """
        le = QLineEdit()

        le.setText(default_text)
        le.setValidator(QDoubleValidator())
        #le.setMaxLength(max_length)
        le.setAlignment(Qt.AlignRight)
        le.setFixedWidth(200)

        return le
    
    def plot_graph(self, k, i):
        ax = self.figures[i].add_subplot(111)
        ax.clear()

        if i == 0:
            # Vanilla
            if self.mode_cp == 0:
                # Call
                ax.plot([0, k], [0, 0],"g")
                x = np.arange(k,k+50,5)
                y = np.arange(10)
                ax.plot(x, y,"g")
            else:
                # Put
                x = np.arange(k+1)
                y = [k-i for i in range(int(k)+1)]
                ax.plot([k, k+k], [0, 0],"g")
                ax.plot(x, y,"g")
        elif i == 1:
            # Cash
            # Calculation mode is Call, else Put 
            ax.plot([1, k, k, k+k] if self.mode_cp == 0 else [k+k,k, k, 1], [0, 0, 1, 1])
        else:
            # Asset
            x = np.arange(k)
            y = np.arange(k, k+k, 1)
            z = np.arange(k, k+k, 1)

            ax.plot(z,y,"w" if self.mode_cp == 1 else "b")
            ax.plot(x, "b" if self.mode_cp == 1 else "w")
        
        # The selected model is not the BS or MCS
        # The graph in question (i) is the cash or asset
        ax.set_visible(not (self.model not in [0,2] and i in [1,2]))
        
        self.graphs[i].draw()

    def get_data(self):
        """Function to get the data
        """

        self.lb_update_progress.setText('Updating')
        self.btn_get_ticker_data.setEnabled(False)
        ticker = self.tb_ticker.text().upper()
        self.tb_ticker.setText(ticker)

        try:
            # Try to get the ticker data and update the fields
            data = get_ticker_data(ticker)
            
            # 'S','K', 'r', 'sigma', 'T', 'q'
            self.tb_v[0].setText(str(data['S']))
            self.tb_v[2].setText(str(data['r']))
            self.tb_v[3].setText(str(data['sigma']))
            self.tb_v[5].setText(str(data['q']))
            self.lb_update_progress.setText('')
        except:
            # Caught when the ticker cannot be found
            self.lb_update_progress.setText('Error')

        self.btn_get_ticker_data.setEnabled(True)
    
    def reset_data(self):
        self.tb_ticker.setText('')
        
        for t in self.tb_v:
            t.setText(str(0))
        self.tb_v[2].setText(str(0.017))
        self.tb_v[4].setText(str(183/365))
    
    def calculate(self):
        """Perform option pricing and binary option calculations
        """

        print('Calculate!')
        for i,v in enumerate(self.tb_v):
            option_values[i] = float(v.text())
        # option_values = [float(v.text()) for i,v in enumerate(self.tb_v)]

        op = BinaryOption()
        op.set_variables(option_values)
        self.outputs = {
            0: op.black_scholes(),
            1: op.binomial_tree(),
            2: op.monte_carlo(get_z(10_000_000)),
            3: op.implicit(),
            4: op.explicit(),
            5: op.crank_nicolson(),
        }

        # Set to calculate cash (mode 1)
        op.set_bin_asset_type(1)
        self.bin1 = {
            0: op.binary_black_scholes(),
            1: (0,0),
            2: op.binary_monte_carlo(SIMULATIONS),
            3: (0,0),
            4: (0,0),
            5: (0,0),
        }

        # Set to calculate asset (mode 2)
        op.set_bin_asset_type(2)
        self.bin2 = {
            0: op.binary_black_scholes(),
            1: (0,0),
            2: op.binary_monte_carlo(SIMULATIONS),
            3: (0,0),
            4: (0,0),
            5: (0,0),
        }

        # Display the result
        self.show_result()
        self.cb_models.currentIndexChanged.connect(self.show_result)
    
    def toggle_cp(self):
        """Toggles the Call/Put button
        """
        self.mode_cp = 1 if self.mode_cp == 0 else 0
        self.btn_cp.setText(OPTION_OPTIONS[self.mode_cp])
        self.show_result()

    def show_result(self):
        """Display the result
        """
        self.model = self.cb_models.currentIndex()
        self.lb_figure[0].setText('${:,.2f}'.format(self.outputs[self.model][self.mode_cp]))
        self.lb_figure[1].setText('${:,.2f}'.format(self.bin1[self.model][self.mode_cp]))
        self.lb_figure[2].setText('${:,.2f}'.format(self.bin2[self.model][self.mode_cp]))
        # print(option_values)
        k = option_values[1]

        # Plot the graphs
        self.plot_graph(k,0) # Vanilla
        self.plot_graph(k,1) # Cash
        self.plot_graph(k,2) # Asset

if __name__=='__main__':
    app = QApplication(sys.argv)
    main = Main()
    main.show()
    sys.exit(app.exec_())

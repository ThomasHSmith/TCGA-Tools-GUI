import sys, os
from PyQt4 import QtCore, QtGui, uic
import pandas as pd
from PyQt4.uic import loadUiType
 
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
import numpy as np
from scipy import stats
import seaborn as sns

qtCreatorFile = "./TCGA_Tools_UI.ui"
 
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
 
class MyApp(QtGui.QMainWindow, Ui_MainWindow):
    global choices_dict
    global data_dir
    global df, df_targets, df_targets_log2, df_targets_z, df_targets_log2_z
    global canvasFull, log2opt, ZScoreOpt, ZCutOffOpt
    ZScoreOpt = False
    log2opt = True
    canvasFull = False
    ZCutOffOpt = False
    
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.find_directory_button.clicked.connect(self.SpecifyRootDataDir)
        self.load_project_button.clicked.connect(self.LoadProjectData)
        self.load_targets_list_button.clicked.connect(self.LoadTargetsListFile)
        self.extract_targets_button.clicked.connect(self.ExtractTargets)
        self.plot_hist_button.clicked.connect(self.PlotHist)
        self.log2_checkbox.stateChanged.connect(self.ToggleLog2)
        self.z_transform_checkbox.stateChanged.connect(self.ToggleZScore)
        self.z_cutoff_checkbox.stateChanged.connect(self.ToggleZCutOff)
        self.generate_heatmap_button.clicked.connect(self.GenerateHeatMap)



        
    def SpecifyRootDataDir(self):
        global choices_dict
        global data_dir
        data_dir = QtGui.QFileDialog.getExistingDirectory(None, "Select Directory", '/Users/TS_MBP', QtGui.QFileDialog.ShowDirsOnly)
        L = os.listdir(data_dir)
        choices_dict = {}
        for item in L:
            choices_dict[item[:4]] =item
        print L
        for item in choices_dict.keys():
            self.project_choice_menu.addItem(item)

    def LoadProjectData(self):
        global df
        self.metadata_display_box.clear()
        PROJECT = str(self.project_choice_menu.currentText())
        r = str(data_dir)
        print r
        self.metadata_display_box.addItem("Loading %s dataset..." % PROJECT)
        PICKLE_PATH = os.path.join(r, choices_dict[PROJECT])
        print choices_dict[PROJECT]
        df = pd.read_pickle(PICKLE_PATH)
 	# Convert TCGA classifier codes
	# 11 - Normal tissue -> -2
	# 01 - Solid tumor -> 0
	# 02 - Recurrent solid tumor -> 1
	# 06 - Mets -> 2
        vals_dict = {'11':'NormalControl', '01':'SolidTumor', '02':'RecurrentSolidTumor', '06':'Metastatic'}
        df['SampleTypeText'] = df.SampleType.apply(lambda x: vals_dict[x])
        for item in df.SampleTypeText.unique():
            msg='%s: %d' % (item, sum(df.SampleTypeText==item))
            self.metadata_display_box.addItem(msg)
        msg='Total: %d\n' % len(df.SampleTypeText)
        self.metadata_display_box.addItem(msg)

        
    def LoadTargetsListFile(self):
        TARGETS_INFILE = str(QtGui.QFileDialog.getOpenFileName(self, 'OpenFile'))
        #targets_dict = {}
        f = open(TARGETS_INFILE, 'r')
        print 'Reading gene targets from: %s' % TARGETS_INFILE
        for line in f:
            words = line.split()
            key =' '.join(words[:-1])
            gene_symbol=words[0]
            ID = words[-1]
            self.gene_targets_editable_list.appendPlainText(gene_symbol + ' ' + ID)
            #targets_dict[key] = ID
        f.close()

    def ExtractTargets(self):
        global df_targets, df_targets_log2, df_targets_log2_z, df_targets_z
        targets_dict = {}
        targets = str(self.gene_targets_editable_list.toPlainText())
        lines = targets.split('\n')
        for line in lines:
            words = line.split()
            key = words[0]
            ID = words[-1]
            targets_dict[key] = ID

        df_targets = pd.DataFrame()
        for key in targets_dict:
            try:
                df_targets.insert(0, key, df[targets_dict[key]])
            except KeyError:
                print '%s (%s) not found in Dataset' % (targets_dict[key], key)
        df_targets_log2 = df_targets.apply(lambda x: np.log2(x+1))
        df_targets_log2_z = df_targets_log2.select_dtypes(exclude=['object']).apply(stats.zscore)
        df_targets_z = df_targets.select_dtypes(exclude=['object']).apply(stats.zscore)
        self.loaded_targets_list.clear()
        for column in df_targets:
            self.loaded_targets_list.addItem(column)

    def ToggleLog2(self):
        global log2opt
        if self.log2_checkbox.isChecked():
            log2opt = True
            print 'Log2 checked'
        else:
            log2opt = False
            print 'Log2 unchecked'

    def ToggleZScore(self):
        global ZScoreOpt
        if self.z_transform_checkbox.isChecked():
            ZScoreOpt = True
            print 'ZScore checked'
        else:
            ZScoreOpt = False
            print 'ZScore unchecked'

    def ToggleZCutOff(self):
        global ZCutOffOpt
        if self.z_cutoff_checkbox.isChecked():
            ZCutOffOpt = True
        else:
            ZCutOffOpt = False
              
    def addmpl(self, fig):
        global canvasFull
        self.rmmpl()
        self.canvas = FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
        canvasFull = True
        
    def rmmpl(self):
        global canvasFull
        if canvasFull == True:
            self.mplvl.removeWidget(self.canvas)
            self.canvas.close()
            canvasFull = False
        
    def addmpl_2(self, fig):
        self.canvas = FigureCanvas(fig)
        self.mplvl_2.addWidget(self.canvas)
        self.canvas.draw()
        
    def PlotHist(self):
        global df_targets, df_targets_z, df_targets_log2, df_targets_log2_z
        target = str(self.loaded_targets_list.currentItem().text())
        fig1 = Figure()
        ax1f1 = fig1.add_subplot(111)
        if (log2opt & ZScoreOpt):
            if ZCutOffOpt:
                Z_CUTOFF = int(self.z_cutoff.text())
                df_targets_log2_z = df_targets_log2_z[~(df_targets_log2_z > Z_CUTOFF).any(axis=1)]
                df_targets_log2_z = df_targets_log2_z[~(df_targets_log2_z < -Z_CUTOFF).any(axis=1)]
                ax1f1.hist(df_targets_log2_z[target])
            else:
                ax1f1.hist(df_targets_log2_z[target])
        elif (log2opt & (not ZScoreOpt)):
            ax1f1.hist(df_targets_log2[target])
        elif ((not log2opt) & ZScoreOpt):
            ax1f1.hist(df_targets_z[target])
        else:
            ax1f1.hist(df_targets[target])
        self.addmpl(fig1)
        print target
        
    def GenerateHeatMap(self):
        fig2 = Figure()
        ax1f2 = fig2.add_subplot(111)
        sns.heatmap(df_targets_log2_z,yticklabels=False, xticklabels=True, ax=ax1f2)
        self.addmpl_2(fig2)

        
        
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = MyApp()
    window.show()
    sys.exit(app.exec_())

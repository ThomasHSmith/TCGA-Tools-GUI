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
import matplotlib.pyplot as plt

qtCreatorFile = "./TCGA_Tools_UI.ui"
 
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
 
class MyApp(QtGui.QMainWindow, Ui_MainWindow):
    global canvasFull, log2opt, ZScoreOpt, ZCutOffOpt, includeControls
    global data_dir, choices_dict
    # Default values for checkbox bools
    ZScoreOpt = False
    data_dir = '/Users/TS_MBP/TCGA_Pickles'
    log2opt = True
    canvasFull = False
    ZCutOffOpt = False
    includeControls = False

    
    def __init__(self):
        global choices_dict
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
        self.transform_data_button.clicked.connect(self.TransformData)
        self.generate_scatterplot_button.clicked.connect(self.GenerateScatterPlot)
        self.lmplot_button.clicked.connect(self.GeneratelmPlot)
        self.calculate_correlation_button.clicked.connect(self.CalculateCorrelation)
        self.include_controls_checkbox.stateChanged.connect(self.ToggleIncludeControls)

        L = os.listdir(data_dir)
        if len(L) > 0:
    	    self.log_display_box.appendPlainText("Looking for data pickles in default directory: %s" % data_dir)
    	    choices_dict = {}
    	    for item in L:
        	    choices_dict[item[:4]] =item
        	    print item
       	    for item in choices_dict.keys():
                self.project_choice_menu.addItem(item)
        else:
            self.log_display_box.appendPlainText("No data pickles found in default directory. Please specifcy directory")


        
    def SpecifyRootDataDir(self):
        global choices_dict
        global data_dir
        data_dir = QtGui.QFileDialog.getExistingDirectory(None, "Select Directory", '/Users/TS_MBP', QtGui.QFileDialog.ShowDirsOnly)
        L = os.listdir(data_dir)
        self.log_display_box.appendPlainText("Looking for data pickles in %s" % data_dir)
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
        self.log_display_box.appendPlainText("Loading %s dataset..." % PROJECT)
        PICKLE_PATH = os.path.join(r, choices_dict[PROJECT])
        print choices_dict[PROJECT]
        df = pd.read_pickle(PICKLE_PATH)
        self.log_display_box.appendPlainText("Loaded %s!" % PROJECT)
        self.metadata_display_box.addItem("%s metadata:" % PROJECT)
        types_list = list(df.SampleType)
        types_set = set(types_list)
        stages_list = list(df.TumorStage)
        stages_set = set(stages_list)
        for sample_type in types_set:
            msg='%s: %d' % (sample_type, types_list.count(sample_type))
            self.metadata_display_box.addItem(msg)
        msg='Total: %d\n' % len(types_list)
        self.metadata_display_box.addItem(msg)

        for stage in stages_set:
            msg='%s: %d' % (stage, stages_list.count(stage))
            self.metadata_display_box.addItem(msg)
        msg='Total: %d\n' % len(stages_list)
        self.metadata_display_box.addItem(msg)



        
    def LoadTargetsListFile(self):
        TARGETS_INFILE = str(QtGui.QFileDialog.getOpenFileName(self, 'OpenFile'))
        f = open(TARGETS_INFILE, 'r')
        print 'Reading gene targets from: %s' % TARGETS_INFILE
        self.log_display_box.appendPlainText("Reading genes targets from: %s" % TARGETS_INFILE)
        for line in f:
            words = line.split()
            key =' '.join(words[:-1])
            gene_symbol=words[0]
            ID = words[-1]
            self.gene_targets_editable_list.appendPlainText(gene_symbol + ' ' + ID)
        f.close()
    
    def ExtractTargets(self):
        global df_targets, df
        self.log_display_box.appendPlainText("Extracting genes from main dataframe")
        targets_dict = {}
        targets = str(self.gene_targets_editable_list.toPlainText())
        lines = targets.split('\n')
        for line in lines:
            words = line.split()
            gene_name = words[0]
            ensID = words[-1]
            targets_dict[ensID] = gene_name

        df_targets = df.copy()
        missing_keys = []
        for key in targets_dict:
            if key in df.columns:
                #df_targets.insert(0, key, df[targets_dict[key]])
                self.log_display_box.appendPlainText("Found %s (%s) in Dataset" % (key, targets_dict[key]))
            else:
                missing_keys.append(key)
                print '%s (%s) not found in Dataset' % (key, targets_dict[key])
                self.log_display_box.appendPlainText("%s (%s) not found in Dataset" % (key, targets_dict[key]))
        # Remove missing keys from targets_dict
        for key in missing_keys:
            targets_dict.pop(key)

        METADATA_COLS = ['PtID','TumorStage','SampleType']
        targets = targets_dict.keys()
        targets.extend(METADATA_COLS)
        df_targets = df_targets[targets]
        df_targets = df_targets.rename(columns=targets_dict)
        
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
        self.toolbar = NavigationToolbar(self.canvas,
                                         self.mplwindow, coordinates=True)
        self.mplvl.addWidget(self.toolbar)
        
    def rmmpl(self):
        global canvasFull
        if canvasFull == True:
            self.mplvl.removeWidget(self.canvas)
            self.mplvl.removeWidget(self.toolbar)
            self.canvas.close()
            self.toolbar.close()
            canvasFull = False

    def TransformData(self):
        global df_targets, df_transformed
        #TARGET_GENE = str(self.loaded_targets_list.currentItem().text())
        df_transformed = df_targets.copy()
        value_cols = list(df_transformed.select_dtypes(exclude=['object']))
        
        
        if log2opt:
            df_transformed[value_cols] = df_transformed[value_cols].apply(lambda x: np.log2(x+1))
            
        if ZScoreOpt:
            df_transformed[value_cols] = df_transformed[value_cols].apply(stats.zscore)
            if ZCutOffOpt:
                Z_CUTOFF = int(self.z_cutoff.text())
                df_transformed = df_transformed[~(df_transformed[value_cols] > Z_CUTOFF).any(axis=1)]
                df_transformed = df_transformed[~(df_transformed[value_cols] < -Z_CUTOFF).any(axis=1)]
        # Insert MetaData
        #vals_dict = {'11':'NormalControl', '01':'SolidTumor', '02':'RecurrentSolidTumor', '06':'Metastatic'}
        #df_targets.insert(len(df_targets.columns),'SampleType',meta_col)
        df_transformed = df_transformed.drop(labels=['PtID','TumorStage'], axis=1)
        SampleToInt = {'Solid Tissue Normal':-3, 'Primary Tumor':0, 'Metastatic':0, 'Recurrent Solid Tumor':0}
        df_transformed['SampleType'] = df_transformed['SampleType'].apply(lambda x: SampleToInt[x])
        ControlTumorConv = {-3:'Control', 0:'Tumor'}
        df_transformed['Tissue'] = df_transformed['SampleType'].apply(lambda x: ControlTumorConv[x])
        self.log_display_box.appendPlainText("Finished transforming data")

    def ToggleIncludeControls(self):
        global includeControls
        if self.include_controls_checkbox.isChecked():
            includeControls = True
        else:
            includeControls = False

            
    def CalculateCorrelation(self):
        self.log_display_box.appendPlainText("Calculating correlation coefficients...")
        TARGET_GENE = str(self.loaded_targets_list.currentItem().text())
        df_temp = df_transformed.copy()
        if includeControls == False:
            df_temp = df_temp[df_temp['Tissue'] != 'Control']
            self.log_display_box.appendPlainText("Dropped Control samples")
        df_temp = df_temp.select_dtypes(exclude=['object'])

        target_col = df_temp[TARGET_GENE]
        df_temp = df_temp.drop(labels=[TARGET_GENE], axis=1)
        df_temp.insert(0, TARGET_GENE, target_col)
        self.correlation_output_box.clear()
        x=list(df_temp[TARGET_GENE])
        header='Gene\t PearsonR (p)\t\t SpearmanR (p)'
        self.correlation_output_box.appendPlainText(header)

        for column in df_temp:
            y = list(df_temp[column])
            pearsonR, pearsonP = stats.pearsonr(x, y)
            spearmanR, spearmanP = stats.spearmanr(x, y)
            line = '%s\t%.3f (%.3E)\t%.3f (%.3E)' % (column, pearsonR, pearsonP, spearmanR, spearmanP)
            self.correlation_output_box.appendPlainText(line)

        self.log_display_box.appendPlainText("Completed correlation calculations.")       
        

    def PlotHist(self):
        TARGET_GENE = str(self.loaded_targets_list.currentItem().text())
        fig1 = Figure()
        ax1f1 = fig1.add_subplot(111)
        ax1f1.hist(df_transformed[TARGET_GENE])
        self.addmpl(fig1)
        print TARGET_GENE

    def GenerateScatterPlot(self):
        selected_genes = self.loaded_targets_list.selectedItems()
        if len(selected_genes) == 2:
            gene_x = str(selected_genes[0].text())
            gene_y = str(selected_genes[1].text())
            sns.jointplot(x=gene_x, y=gene_y, data=df_transformed, kind='reg', size=10, space=0)
            plt.show()
        else:
            print 'Select only 2 genes'

    def GeneratelmPlot(self):
        selected_genes = self.loaded_targets_list.selectedItems()
        if len(selected_genes) == 2:
            gene_x = str(selected_genes[0].text())
            gene_y = str(selected_genes[1].text())
            sns.lmplot(x=gene_x, y=gene_y, data=df_transformed, hue='Tissue',palette='Set1', size=10)
            plt.show()
        else:
            print 'Select only 2 genes'
       
    def GenerateHeatMap(self):
        TARGET_GENE = str(self.loaded_targets_list.currentItem().text())
        df_temp = df_transformed.copy()
        target_col = df_temp[TARGET_GENE]
        df_temp = df_temp.drop(labels=[TARGET_GENE], axis=1)
        df_temp.insert(0, TARGET_GENE, target_col)
        df_temp = df_temp.sort_values(TARGET_GENE, ascending=False)

        sns.heatmap(df_temp.select_dtypes(exclude=['object']),yticklabels=False, xticklabels=True)
        plt.show()        
        # For embedded:
        #fig1 = Figure()
        #ax1f1 = fig1.add_subplot(111)
        #sns.heatmap(df_temp.select_dtypes(exclude=['object']),yticklabels=False, xticklabels=True, ax=ax1f1)
        #fig1.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.04)
        #self.addmpl(fig1)

        
        
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = MyApp()
    window.show()
    sys.exit(app.exec_())

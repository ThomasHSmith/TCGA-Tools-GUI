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
    global canvasFull, log2opt, ZScoreOpt, ZCutOffOpt, controlsOnly, tumorsOnly
    global data_dir, choices_dict, dataTransformed, get_sample_choice_df
    # Default values for checkbox bools
    ZScoreOpt = False
    data_dir = '/Users/TS_MBP/TCGA_Pickles'
    log2opt = True
    canvasFull = False
    ZCutOffOpt = False
    controlsOnly = False
    tumorsOnly = False
    dataTransformed = False


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
        self.export_to_excel_button.clicked.connect(self.ExportToExcel)
        self.save_heatmap_button.clicked.connect(self.SaveHeatMap)
        self.dataDirMenuButton.connect(self.SpecifyRootDataDir)
        self.saveToExcelMenuButton.connect(self.ExportToExcel)
        self.exitMenuButton.triggered.connect(self.close)




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
        global df, df_ctrls, df_tumors
        global all_samples_label, ctrl_samples_label, tumor_samples_label
        self.correlation_output_box.clear()
        self.metadata_display_box.clear()
        self.sample_type_menu.clear()
        self.stage_selection_box.clear()
        PROJECT = str(self.project_choice_menu.currentText())
        r = str(data_dir)
        print r
        self.log_display_box.appendPlainText("Loading %s dataset..." % PROJECT)
        PICKLE_PATH = os.path.join(r, choices_dict[PROJECT])
        print choices_dict[PROJECT]
        df = pd.read_pickle(PICKLE_PATH)

        self.log_display_box.appendPlainText("Loaded %s!" % PROJECT)
        self.metadata_display_box.addItem("%s metadata:" % PROJECT)

        sample_types_dict = dict(df.SampleType.value_counts())
        df['TumorStage'] = df.TumorStage.apply(lambda x: 'na' if ( (x=='not reported')| (x=='stage x') ) else x)
        tumor_stages_dict = dict(df.TumorStage.value_counts())
        for key in sample_types_dict:
            msg='%s: %d' % (key, sample_types_dict[key])
            self.metadata_display_box.addItem(msg)
        msg='Total: %d\n' % len(df)
        self.metadata_display_box.addItem(msg)

        for key in tumor_stages_dict:
            msg='%s: %d' % (key, tumor_stages_dict[key])
            self.metadata_display_box.addItem(msg)
        msg='Total: %d\n' % len(df)
        self.metadata_display_box.addItem(msg)

        na_val, stage1_val, stage2_val, stage3_val, stage4_val = range(0,5)
        stages_dict = { 'stage i': stage1_val, 'stage ia':stage1_val, 'stage ib':stage1_val,
                        'stage iia':stage2_val, 'stage iib': stage2_val, 'stage ii': stage2_val,
                        'stage iiia': stage3_val, 'stage iiic':stage3_val, 'stage iiib':stage3_val,
                        'stage iii':stage3_val, 'stage iv':stage4_val, 'na':'na'}
        df['TumorStageStr'] = df['TumorStage'].apply(lambda x: str(stages_dict[x]))
        stages_dict = dict(df.TumorStageStr.value_counts())
        stages_str_dict = dict(df.TumorStageStr.value_counts())

        for key in stages_str_dict:
            msg='Stage %s: %s' % (key, stages_str_dict[key])
            self.stage_selection_box.addItem(msg)


        df_ctrls = df[df.SampleType == 'Solid Tissue Normal'].copy()
        df_tumors = df[df.SampleType != 'Solid Tissue Normal'].copy()

        all_samples_label = 'All samples (%d)' % len(df)
        ctrl_samples_label = 'Controls only (%d)' % len(df_ctrls)
        tumor_samples_label = 'Tumors only (%d)' % len(df_tumors)

        self.sample_type_menu.addItem(all_samples_label)
        self.sample_type_menu.addItem(ctrl_samples_label)
        self.sample_type_menu.addItem(tumor_samples_label)

    def LoadTargetsListFile(self):
        TARGETS_INFILE = str(QtGui.QFileDialog.getOpenFileName(self, 'OpenFile'))
        f = open(TARGETS_INFILE, 'r')
        print 'Reading gene targets from: %s' % TARGETS_INFILE
        self.correlation_output_box.clear()
        self.gene_targets_editable_list.clear()
        self.log_display_box.appendPlainText("Reading genes targets from: %s" % TARGETS_INFILE)
        for line in f:
            words = line.split()
            key =' '.join(words[:-1])
            gene_symbol=words[0]
            ID = words[-1]
            self.gene_targets_editable_list.appendPlainText(gene_symbol + ' ' + ID)
        f.close()

    def ExtractTargets(self):
        global df_targets, df, df_targets_controls, df_targets_tumors
        self.correlation_output_box.clear()
        self.log_display_box.appendPlainText("Extracting genes from main dataframe")
        targets_dict = {}
        IDs_ordered = []
        targets = str(self.gene_targets_editable_list.toPlainText())
        lines = targets.split('\n')
        for line in lines:
            words = line.split()
            gene_name = ' '.join(words[:1]).strip().rstrip()
            ensID = words[-1]
            targets_dict[ensID] = gene_name
            IDs_ordered.append(ensID)

        METADATA_COLS = ['PtID','TumorStage','TumorStageStr','SampleType']
        df_targets = df[IDs_ordered+METADATA_COLS].copy()
        df_targets.rename(columns=targets_dict, inplace=True)

        df_targets_controls = df_ctrls[IDs_ordered+METADATA_COLS].copy()
        df_targets_controls.rename(columns=targets_dict, inplace=True)

        df_targets_tumors = df_tumors[IDs_ordered+METADATA_COLS].copy()
        df_targets_tumors.rename(columns=targets_dict, inplace=True)

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
        global df_targets, df_transformed, dataTransformed
        global df_targets_controls, df_targets_tumors
        global df_controls_trans, df_tumors_trans
        #TARGET_GENE = str(self.loaded_targets_list.currentItem().text())
        df_transformed = df_targets.copy()
        df_controls_trans = df_targets_controls.copy()
        df_tumors_trans = df_targets_tumors.copy()

        value_cols = list(df_transformed.select_dtypes(exclude=['object']))

        if log2opt:
            df_transformed[value_cols] = df_transformed[value_cols].apply(lambda x: np.log2(x+1))
            df_controls_trans[value_cols] = df_controls_trans[value_cols].apply(lambda x: np.log2(x+1))
            df_tumors_trans[value_cols] = df_tumors_trans[value_cols].apply(lambda x: np.log2(x+1))
            dataTransformed = True

        if ZScoreOpt:
            dataTransformed = True
            df_transformed[value_cols] = df_transformed[value_cols].apply(stats.zscore)
            df_controls_trans[value_cols] = df_controls_trans[value_cols].apply(stats.zscore)
            df_tumors_trans[value_cols] = df_tumors_trans[value_cols].apply(stats.zscore)
            if ZCutOffOpt:
                Z_CUTOFF = int(self.z_cutoff.text())
                df_transformed = df_transformed[~(df_transformed[value_cols] > Z_CUTOFF).any(axis=1)]
                df_transformed = df_transformed[~(df_transformed[value_cols] < -Z_CUTOFF).any(axis=1)]

                df_controls_trans = df_controls_trans[~(df_controls_trans[value_cols] > Z_CUTOFF).any(axis=1)]
                df_controls_trans = df_controls_trans[~(df_controls_trans[value_cols] < -Z_CUTOFF).any(axis=1)]

                df_tumors_trans = df_tumors_trans[~(df_tumors_trans[value_cols] > Z_CUTOFF).any(axis=1)]
                df_tumors_trans = df_tumors_trans[~(df_tumors_trans[value_cols] < -Z_CUTOFF).any(axis=1)]

        # Insert MetaData
        #vals_dict = {'11':'NormalControl', '01':'SolidTumor', '02':'RecurrentSolidTumor', '06':'Metastatic'}
        #df_targets.insert(len(df_targets.columns),'SampleType',meta_col)

        df_transformed = df_transformed.drop('PtID', axis=1)
        df_controls_trans = df_controls_trans.drop('PtID', axis=1)
        df_tumors_trans = df_tumors_trans.drop('PtID', axis=1)

        SampleToInt = {'Solid Tissue Normal':-3, 'Primary Tumor':0, 'Metastatic':0, 'Recurrent Solid Tumor':0, 'Recurrent Tumor':0}
        df_transformed['SampleType'] = df_transformed['SampleType'].apply(lambda x: SampleToInt[x])
        df_controls_trans['SampleType'] = df_controls_trans['SampleType'].apply(lambda x: SampleToInt[x])
        df_tumors_trans['SampleType'] = df_tumors_trans['SampleType'].apply(lambda x: SampleToInt[x])

        ControlTumorConv = {-3:'Control', 0:'Tumor'}
        df_transformed['Tissue'] = df_transformed['SampleType'].apply(lambda x: ControlTumorConv[x])
        df_controls_trans['Tissue'] = df_controls_trans['SampleType'].apply(lambda x: ControlTumorConv[x])
        df_tumors_trans['Tissue'] = df_tumors_trans['SampleType'].apply(lambda x: ControlTumorConv[x])

        self.log_display_box.appendPlainText("Finished transforming data")

    def get_sample_choice_df(self):
        stage_choices = [(str(x.text()).split(':')[0].split(' ')[1]) for x in self.stage_selection_box.selectedItems()]
        print stage_choices
        if len(stage_choices) == 0:
            stage_choices = ['na', '1', '2', '3', '4']

        sample_choice = str(self.sample_type_menu.currentText())
        if ((sample_choice == '') | (sample_choice == all_samples_label)):
            if dataTransformed:
                return ( df_transformed[df_transformed.TumorStageStr.isin(stage_choices)].copy(), 'all samples (transformed')
            else:
                return (df_targets[df_targets.TumorStageStr.isin(stage_choices)].copy(), 'all samples (untransformed)')
        elif sample_choice == ctrl_samples_label:
            if dataTransformed:
                return (df_controls_trans[df_controls_trans.TumorStageStr.isin(stage_choices)].copy(), 'control samples only (transformed)')
            else:
                return (df_targets_controls[df_targets_controls.TumorStageStr.isin(stage_choices)].copy(), 'control samples only (untransformed)')
        elif sample_choice == tumor_samples_label:
            if dataTransformed:
                return (df_tumors_trans[df_tumors_trans.TumorStageStr.isin(stage_choices)].copy(), 'tumor samples only (transformed)')
            else:
                return (df_targets_tumors[df_targets_tumors.TumorStageStr.isin(stage_choices)].copy(), 'tumor samples only (untransformed)')



    def CalculateCorrelation(self):
        TARGET_GENE = str(self.loaded_targets_list.currentItem().text())
        df_temp, samp_choice = get_sample_choice_df(self)
        self.correlation_output_box.clear()

        data_summary_header = 'Sample Types Included\tNumber of Samples'
        self.correlation_output_box.appendPlainText(data_summary_header)
        sample_types_dict = dict(df_temp.SampleType.value_counts())
        cum_sum = 0
        for key in sample_types_dict:
            msg='\t%s\t\t%d' % (key, sample_types_dict[key])
            cum_sum = cum_sum + int(sample_types_dict[key])
            self.correlation_output_box.appendPlainText(msg)
        msg='Total:\t\t%d\n' % cum_sum
        self.correlation_output_box.appendPlainText(msg)

        data_summary_header = 'Tumor Stages Included\tNumber of Samples'
        self.correlation_output_box.appendPlainText(data_summary_header)
        tumor_stages_dict = dict(df_temp.TumorStageStr.value_counts())
        cum_sum = 0
        for key in tumor_stages_dict:
            msg='\tStage %s\t\t%d' % (key, tumor_stages_dict[key])
            cum_sum = cum_sum + int(tumor_stages_dict[key])
            self.correlation_output_box.appendPlainText(msg)
        msg='Total:\t\t%d\n' % cum_sum
        self.correlation_output_box.appendPlainText(msg)
        msg = '-'*30+'\n'
        self.correlation_output_box.appendPlainText(msg)

        TARGET_GENE = str(self.loaded_targets_list.currentItem().text())
        df_temp = df_temp.select_dtypes(exclude=['object'])
        target_col = df_temp[TARGET_GENE]
        df_temp = df_temp.drop(labels=[TARGET_GENE], axis=1)
        df_temp.insert(0, TARGET_GENE, target_col)
        x=list(df_temp[TARGET_GENE])
        header='Gene\t PearsonR (p)\t\t SpearmanR (p)'
        self.correlation_output_box.appendPlainText(header)

        for column in df_temp:
            y = list(df_temp[column])
            pearsonR, pearsonP = stats.pearsonr(x, y)
            spearmanR, spearmanP = stats.spearmanr(x, y)
            line = '%s\t%.3f (%.3E)\t%.3f (%.3E)' % (column, pearsonR, pearsonP, spearmanR, spearmanP)
            self.correlation_output_box.appendPlainText(line)
        self.log_display_box.appendPlainText(("Calculated correlation coefficients for %s from %s..." % (TARGET_GENE, samp_choice)) )


    def PlotHist(self):
        TARGET_GENE = str(self.loaded_targets_list.currentItem().text())
        fig1 = Figure()
        ax1f1 = fig1.add_subplot(111)
        df_temp, samp_choice = get_sample_choice_df(self)
        print 'Loaded %s' % samp_choice
        ax1f1.hist(df_temp[TARGET_GENE])
        self.addmpl(fig1)
        print TARGET_GENE

    def GenerateScatterPlot(self):
        selected_genes = self.loaded_targets_list.selectedItems()
        if len(selected_genes) == 2:
            gene_x = str(selected_genes[0].text())
            gene_y = str(selected_genes[1].text())
            df_temp, samp_choice = get_sample_choice_df(self)
            self.log_display_box.appendPlainText(("Generated scatter plot for %s and %s from %s" % (gene_x, gene_y, samp_choice)) )
            sns.jointplot(x=gene_x, y=gene_y, data=df_temp, kind='reg', size=10, space=0)
            plt.show()

#            if dataTransformed:
#                sns.jointplot(x=gene_x, y=gene_y, data=df_transformed, kind='reg', size=10, space=0)
#            else:
#                sns.jointplot(x=gene_x, y=gene_y, data=df_targets, kind='reg', size=10, space=0)
#            plt.show()
        else:
            print 'Select only 2 genes'
            self.log_display_box.appendPlainText('Select only 2 genes')

    def GeneratelmPlot(self):
        selected_genes = self.loaded_targets_list.selectedItems()
        if len(selected_genes) == 2:
            gene_x = str(selected_genes[0].text())
            gene_y = str(selected_genes[1].text())
            df_temp, samp_choice = get_sample_choice_df(self)
            self.log_display_box.appendPlainText(("Generated lm plot for %s and %s from %s" % (gene_x, gene_y, samp_choice)) )
            sns.lmplot(x=gene_x, y=gene_y, data=df_transformed, hue='TumorStage',palette='Set1', size=10)
            plt.show()

#            if dataTransformed:
#                sns.lmplot(x=gene_x, y=gene_y, data=df_transformed, hue='TumorStage',palette='Set1', size=10)
#            else:
#                sns.lmplot(x=gene_x, y=gene_y, data=df_targets, hue='TumorStage',palette='Set1', size=10)
#            plt.show()
        else:
            print 'Select only 2 genes'
            self.log_display_box.appendPlainText('Select only 2 genes')


    def GenerateHeatMap(self):
        TARGET_GENE = str(self.loaded_targets_list.currentItem().text())
        df_temp, samp_choice = get_sample_choice_df(self)
#        if dataTransformed:
#            df_temp = df_transformed.copy()
#        else:
#            df_temp = df_targets.copy()
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


    def SaveHeatMap(self):
        fileName = str(QtGui.QFileDialog.getSaveFileName(self, 'Save as...', '/Users/TS_MBP/Documents/GitHub/TCGA-Tools-GUI/', selectedFilter='*.svg'))

        TARGET_GENE = str(self.loaded_targets_list.currentItem().text())
        df_temp, samp_choice = get_sample_choice_df(self)

#        if dataTransformed:
#            df_temp = df_transformed.copy()
#        else:
#            df_temp = df_targets.copy()
        target_col = df_temp[TARGET_GENE]
        df_temp = df_temp.drop(labels=[TARGET_GENE], axis=1)
        df_temp.insert(0, TARGET_GENE, target_col)
        df_temp = df_temp.sort_values(TARGET_GENE, ascending=False)


        #fig_height = int( (len(df_temp))/5)
        #fig = plt.figure(figsize=(3,fig_height))
        #fig = Figure(figsize=(3, fig_height), dpi=300)
        #fig1 = Figure(figsize=(3,fig_height))
        ax = sns.heatmap(df_temp.select_dtypes(exclude=['object']), yticklabels=False, xticklabels=True)
        file_format = fileName.split('.')[-1]
        fig = ax.get_figure()
        fig.savefig(fileName, format=file_format, transparent=True, bbox_inchest='tight')
        #hm.savefig(fileName, dpi=300, format=file_format, transparent=True, bbox_inches='tight')
        self.log_display_box.appendPlainText("Saved heatmap as %s" % fileName)
#        baseFileName = fileName.split('.')[0]
#        excelFileName = '%s.xlsx' % baseFileName
#        df_temp.to_excel(excelFileName)
        plt.close(fig)

    def ExportToExcel(self):
        fileName = str(QtGui.QFileDialog.getSaveFileName(self, 'Save as...', '/Users/TS_MBP/Documents/GitHub/TCGA-Tools-GUI/', selectedFilter='*.xlsx'))
        TARGET_GENE = str(self.loaded_targets_list.currentItem().text())
        df_temp, samp_choice = get_sample_choice_df(self)

#        if dataTransformed:
#            df_temp = df_transformed.copy()
#        else:
#            df_temp = df_targets.copy()
        target_col = df_temp[TARGET_GENE]
        df_temp = df_temp.drop(labels=[TARGET_GENE], axis=1)
        df_temp.insert(0, TARGET_GENE, target_col)
        df_temp = df_temp.sort_values(TARGET_GENE, ascending=False)
        baseFileName = fileName.split('.')[0]
        excelFileName = '%s.xlsx' % baseFileName
        df_temp.to_excel(excelFileName)
        self.log_display_box.appendPlainText(("Saved %s data to Excel table as %s" % (samp_choice, excelFileName)))

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = MyApp()
    window.show()
    sys.exit(app.exec_())

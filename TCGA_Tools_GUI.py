import sys, os, webbrowser, bs4, requests
from PyQt4 import QtCore, QtGui, uic
from PyQt4.QtGui import *
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
import re
import urllib2
from urllib2 import HTTPError, URLError
import xml.dom.minidom
from urllib2 import HTTPError
from operator import itemgetter

qtCreatorFile = "./UI_layout.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class MyApp(QtGui.QMainWindow, Ui_MainWindow):
    global canvasFull, log2Opt, ZScoreOpt, ZCutOffOpt, controlsOnly, tumorsOnly
    global data_dir, choices_dict, dataTransformed, get_sample_choice_df
    global haveScatter, haveHeatmap, haveHistogram

    # Default values for checkbox bools
    ZScoreOpt = False
    data_dir = '/Users/TS_MBP/TCGA_Pickles'
    log2Opt = True
    canvasFull = False
    ZCutOffOpt = False
    controlsOnly = False
    tumorsOnly = False
    dataTransformed = False
    haveScatter, haveHeatmap, haveHistogram = False, False, False

    def __init__(self):
        global choices_dict
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        w = QWidget()
        self.setupUi(self)

        self.find_directory_button.clicked.connect(self.SpecifyRootDataDir)
        self.load_project_button.clicked.connect(self.LoadProjectData)
        self.load_targets_list_button.clicked.connect(self.LoadTargetsListFile)
        self.extract_targets_button.clicked.connect(self.ExtractTargets)
        self.plot_hist_button.clicked.connect(self.PlotHist)
        self.log2_button.clicked.connect(self.Log2Transform)
        self.z_score_button.clicked.connect(self.ZScore)
        self.z_cutoff_box.valueChanged.connect(self.ZCutOff)
        self.apply_button.clicked.connect(self.apply_clicked)
        self.generate_heatmap_button.clicked.connect(self.GenerateHeatMap)
        self.reset_data_button.clicked.connect(self.ResetData)
        self.generate_scatterplot_button.clicked.connect(self.GenerateScatterPlot)
        self.calculate_correlation_button.clicked.connect(self.CalculateCorrelation)
        self.export_to_excel_button.clicked.connect(self.ExportToExcel)
        self.quitButton.clicked.connect(self.quitApp)
        self.gene_search_button.clicked.connect(self.SearchForGeneIDs)
        self.add_to_list_button.clicked.connect(self.AddResultToGeneList)
        self.clear_results_button.clicked.connect(self.ClearSearchResults)
        self.open_ensembl_page.clicked.connect(self.OpenEnsemblPage)
        self.open_uniprot_page.clicked.connect(self.OpenUniprotPage)
        self.select_all_button.clicked.connect(self.SelectAllResults)

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

    def quitApp(self):
        sys.exit()

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
        global df, df_ctrls, df_tumors, ZScoreOpt, log2Opt, dataTransformed, haveData
        global all_samples_label, ctrl_samples_label, tumor_samples_label

        self.loaded_targets_list.clear()
        self.data_status_box.clear()
        haveData, dataTransformed = False, False
        log2Opt, ZScoreOpt = False, False

        self.data_status_box.addItem('Log2 transformed:\tNo')
        self.data_status_box.addItem('Z-score:\tNo')
        self.data_status_box.addItem('Z-score cut-off:\tNo')

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
                        'stage iii':stage3_val, 'stage iv':stage4_val, 'stage iva':stage4_val,
                        'stage ivb': stage4_val, 'stage ivc':stage4_val, 'na':'na'}
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

    def OpenEnsemblPage(self, ensID):
        if self.search_results_list.count() == 0:
            return
        for item in self.search_results_list.selectedItems():
            result = str(item.text())
            print result
            words = result.split()
            gene_symbol = words[:-1]
            ID = words[-1]
            url = 'http://www.ensembl.org/id/%s' % ID
            new = 2
            webbrowser.open(url,new=new)

    def OpenUniprotPage(self, ensID):
        if self.search_results_list.count() == 0:
            return
        for item in self.search_results_list.selectedItems():
            result = str(item.text())
            print result
            words = result.split()
            gene_symbol = words[:-1]
            ID = words[-1]
            uniprot_IDs = self.get_uniprot_IDs(ID)
            if len(uniprot_IDs) > 0:
                for ID in uniprot_IDs:
                    url = 'http://www.uniprot.org/uniprot/%s' % ID
                    new = 2
                    webbrowser.open(url,new=new)
            else: show_msg_box('UniProt page not found','Could not find UniProt IDs associated with the given Ensembl ID.  Double check EnsemblID.')


    def ExtractTargets(self):
        global df_targets, df, df_targets_controls, df_targets_tumors, haveData
        if haveData == True:
            self.ResetData()
        self.correlation_output_box.clear()
        self.log_display_box.appendPlainText("Extracting genes from main dataframe")
        targets_dict = {}
        IDs_ordered = []
        not_found_dict = {}
        targets = str(self.gene_targets_editable_list.toPlainText())
        lines = targets.split('\n')
        for line in lines:
            words = line.split()
            gene_name = ' '.join(words[:1]).strip().rstrip()
            ensID = words[-1]
            if ensID in df.columns:
                targets_dict[ensID] = gene_name
                IDs_ordered.append(ensID)
            else:
                not_found_dict[gene_name] = ensID
        haveData = True
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

        if len(not_found_dict) > 0:
            missing = ''
            for key in not_found_dict:
                missing = missing+key+': '+not_found_dict[key]+'\n'
            self.show_msg_box('Some Ensembl IDs not found in dataset', missing)

    def SearchForGeneIDs(self):
        self.search_results_list.clear()
        search_targets = str(self.gene_search_input.toPlainText())
        items = search_targets.split(',')
        targets = [x.strip().rstrip() for x in items]
        for target in targets:
            print target
            res_dict = self.get_ensIDs(target)
            if len(res_dict) > 0:
                for key in res_dict:
                    msg = '%s %s' % (key, res_dict[key])
                    self.search_results_list.addItem(msg)

    def SelectAllResults(self):
        if self.search_results_list.count() == 0:
            return
        for item in xrange(self.search_results_list.count()):
            self.search_results_list.setCurrentItem(item)

    def AddResultToGeneList(self):
        for item in self.search_results_list.selectedItems():
            result = str(item.text())
        #results_to_add = [ str(x.text()) for x in self.search_results_list.selectedItems()]
        #for result in results_to_add:
            if len(result) > 1:
                words = result.split()
                name=' '.join(words[:-1])
                ID=words[-1]
                curr_targets = str(self.gene_targets_editable_list.toPlainText())
                lines = curr_targets.split('\n')
                IDs_ordered = []
                if len(lines) > 1:
                    for line in lines:
                        curr_words = line.split()
                        ensID = curr_words[-1]
                        IDs_ordered.append(ensID)
                    if len(IDs_ordered) > 0:
                        if ID in IDs_ordered:
                            self.search_results_list.takeItem(self.search_results_list.row(item))
                        else:
                            self.gene_targets_editable_list.appendPlainText(name + ' ' + ID)
                            self.search_results_list.takeItem(self.search_results_list.row(item))
                else:
                    self.gene_targets_editable_list.appendPlainText(name + ' ' + ID)
                    self.search_results_list.takeItem(self.search_results_list.row(item))
            else:
                print 'Invalid selection'

    def ClearSearchResults(self):
        self.search_results_list.clear()

    def get_uniprot_IDs(self, input_name):
        url = 'http://www.ensembl.org/id/%s' % input_name
        page = requests.get(url)
        soup = bs4.BeautifulSoup(page.text)
        a = soup.findAll('td')
        IDs = []
        for line in a:
            words = str(line).split('http://www.uniprot.org/')
            for word in words:
                if word.find('uniprot/') != -1:
                    IDs.append(word.strip().split('uniprot/')[1].split('"')[0])
        return IDs

    def get_ensIDs(self, input_name):
        try:
            url = 'http://www.uniprot.org/uniprot/%s.xml' % input_name
            fp = urllib2.urlopen(url)
            fp.close()
            IDs = []
            IDs.append(input_name)
        except HTTPError:
            url = 'http://www.uniprot.org/uniprot/?query=gene:%s+AND+organism:9606&format=tab&columns=id' % input_name
            data = urllib2.urlopen(url)
            IDs = [x.rstrip() for x in data.readlines()[1:]]
        ids_dict = {}
        for uniprot_ID in IDs:
            url = 'http://www.uniprot.org/uniprot/%s.xml' % uniprot_ID
            fp = urllib2.urlopen(url)
            doc = xml.dom.minidom.parse(fp)
            fp.close()
            ensIDs = []
            dbrefs = doc.getElementsByTagName('dbReference')
            names_tag = doc.getElementsByTagName('gene')[0]
            gene_name = str(names_tag.getElementsByTagName('name')[0].childNodes[0].data)
            for ref in dbrefs:
                val = str(ref.getAttributeNode('id').value)
                if len(re.findall('ENSG\d.{1,15}', val)) > 0:
                    ensIDs.append(val)
            ensIDs = set(ensIDs)
            ensIDs = list(ensIDs)
            if len(ensIDs) > 1:
                print '%d unique ENSG IDs found for %s (%s)' % (len(ensIDs), gene_name, uniprot_ID)
            if len(ensIDs) == 1:
                ids_dict[gene_name] = ensIDs[0]
        return ids_dict


    def Log2Transform(self):
        global log2Opt
        log2Opt = True
        self.TransformData()

    def ZScore(self):
        global ZScoreOpt
        ZScoreOpt = True
        self.TransformData()

    def ZCutOff(self):
        global ZCutOffOpt
        Z_CUTOFF = float(self.z_cutoff_box.value())
        if Z_CUTOFF != 0.0:
            ZCutOffOpt = True
            self.TransformData()
        else:
            ZCutOffOpt = False
            self.TransformData()

    def addmpl(self, fig):
        global canvasFull
        if canvasFull:
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

    def show_msg_box(self, msg, details):
        mbox = QMessageBox(w)
        mbox.setIcon(QMessageBox.Warning)
        mbox.setText(msg)
        mbox.setDetailedText(details)
#        mbox.setInformativeText(details)
        mbox.setStandardButtons(QMessageBox.Ok)
        mbox.exec_()

    def ResetData(self):
        global log2Opt, ZScoreOpt, ZCutOffOpt, dataTransformed
        dataTransformed, log2Opt, ZScoreOpt, ZCutOffOpt = False, False, False, False
        log2Msg = 'Log2 transformed:\tNo'
        ZScoreMsg = 'Z-score:\tNo'
        ZCutOffMsg = 'Z-score cut-off:\tNo'
        self.data_status_box.clear()
        self.data_status_box.addItem(log2Msg)
        self.data_status_box.addItem(ZScoreMsg)
        self.data_status_box.addItem(ZCutOffMsg)
        self.z_cutoff_box.setValue(0.0)
        self.log_display_box.appendPlainText("Raw data restored")
        self.TransformData()

    def UpdateStagesList(self, df):
        stages_str_dict = dict(df.TumorStageStr.value_counts())
        self.stage_selection_box.clear()
        for key in stages_str_dict:
            msg='Stage %s: %s' % (key, stages_str_dict[key])
            self.stage_selection_box.addItem(msg)

    def apply_clicked(self):
        if haveScatter:
            self.GenerateScatterPlot()
            return
        elif haveHeatmap:
            self.GenerateHeatMap()
            return
        elif haveHistogram:
            self.PlotHist()

    def TransformData(self):
        global df_targets, df_transformed, dataTransformed
        global df_targets_controls, df_targets_tumors
        global df_controls_trans, df_tumors_trans
        global log2Opt, ZScoreOpt, ZCutOffOpt, dataTransformed

        df_transformed = df_targets.copy()
        df_controls_trans = df_targets_controls.copy()
        df_tumors_trans = df_targets_tumors.copy()

        value_cols = list(df_transformed.select_dtypes(exclude=['object']))

        if log2Opt == True:
            df_transformed[value_cols] = df_transformed[value_cols].apply(lambda x: np.log2(x+1))
            df_controls_trans[value_cols] = df_controls_trans[value_cols].apply(lambda x: np.log2(x+1))
            df_tumors_trans[value_cols] = df_tumors_trans[value_cols].apply(lambda x: np.log2(x+1))
            dataTransformed = True

        if ZScoreOpt == True:
            dataTransformed = True
            df_transformed[value_cols] = df_transformed[value_cols].apply(stats.zscore)
            df_controls_trans[value_cols] = df_controls_trans[value_cols].apply(stats.zscore)
            df_tumors_trans[value_cols] = df_tumors_trans[value_cols].apply(stats.zscore)
            if ZCutOffOpt == True:
                Z_CUTOFF = float(self.z_cutoff_box.value())
                if Z_CUTOFF != 0.0:
                    df_transformed = df_transformed[~(df_transformed[value_cols] > Z_CUTOFF).any(axis=1)]
                    df_transformed = df_transformed[~(df_transformed[value_cols] < -Z_CUTOFF).any(axis=1)]

                    df_controls_trans = df_controls_trans[~(df_controls_trans[value_cols] > Z_CUTOFF).any(axis=1)]
                    df_controls_trans = df_controls_trans[~(df_controls_trans[value_cols] < -Z_CUTOFF).any(axis=1)]

                    df_tumors_trans = df_tumors_trans[~(df_tumors_trans[value_cols] > Z_CUTOFF).any(axis=1)]
                    df_tumors_trans = df_tumors_trans[~(df_tumors_trans[value_cols] < -Z_CUTOFF).any(axis=1)]

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

        if log2Opt == True:
            log2Msg = 'Log2 transformed:\tYes'
        elif log2Opt == False:
            log2Msg = 'Log2 transformed:\tNo'

        if ZScoreOpt == True:
            ZScoreMsg = 'Z-score:\tYes'
            if ZCutOffOpt == True:
                Z_CUTOFF = float(self.z_cutoff_box.value())
                if Z_CUTOFF == 0.0:
                    ZCutOffMsg = 'Z-score cut-off:\tNo'
                    ZCutOffOpt = False
                else:
                    ZCutOffMsg = 'Z-score cut-off:\t%f' % Z_CUTOFF
        elif ZScoreOpt == False:
            ZScoreMsg = 'Z-score:\tNo'
            ZCutOffMsg = 'Z-score cut-off:\tNo'
            ZCutOffOpt = False

        self.data_status_box.clear()
        self.data_status_box.addItem(log2Msg)
        self.data_status_box.addItem(ZScoreMsg)
        self.data_status_box.addItem(ZCutOffMsg)


    def get_sample_choice_df(self, force_raw=False):
        global dataTransformed
        stage_choices = [(str(x.text()).split(':')[0].split(' ')[1]) for x in self.stage_selection_box.selectedItems()]
        print stage_choices
        if len(stage_choices) == 0:
            stage_choices = ['na', '1', '2', '3', '4']

        sample_choice = str(self.sample_type_menu.currentText())
        if ((sample_choice == '') | (sample_choice == all_samples_label)):
            if dataTransformed == True:
                if force_raw == False:
                    df_out = df_transformed[df_transformed.TumorStageStr.isin(stage_choices)].copy()
                    #self.UpdateStagesList(df_out)
                    return (df_out, 'all samples (transformed)')
                elif force_raw == True:
                    df_out = df_targets[df_targets.TumorStageStr.isin(stage_choices)].copy()
                    #self.UpdateStagesList(df_out)
                    return (df_out, 'all samples (untransformed)')
            elif dataTransformed == False:
                df_out = df_targets[df_targets.TumorStageStr.isin(stage_choices)].copy()
                #self.UpdateStagesList(df_out)
                return (df_out, 'all samples (untransformed)')

        elif sample_choice == ctrl_samples_label:
            if dataTransformed == True:
                if force_raw == False:
                    df_out = df_controls_trans[df_controls_trans.TumorStageStr.isin(stage_choices)].copy()
                    #self.UpdateStagesList(df_out)
                    return (df_controls_trans[df_controls_trans.TumorStageStr.isin(stage_choices)].copy(), 'control samples only (transformed)')
                elif force_raw == True:
                    df_out = df_targets_controls[df_targets_controls.TumorStageStr.isin(stage_choices)].copy()
                    #self.UpdateStagesList(df_out)
                    return (df_out, 'control samples only (untransformed)')
            elif dataTransformed == False:
                df_out = df_targets_controls[df_targets_controls.TumorStageStr.isin(stage_choices)].copy()
                #self.UpdateStagesList(df_out)
                return (df_out, 'control samples only (untransformed)')

        elif sample_choice == tumor_samples_label:
            if dataTransformed == True:
                if force_raw == False:
                    df_out = df_tumors_trans[df_tumors_trans.TumorStageStr.isin(stage_choices)].copy()
                    #self.UpdateStagesList(df_out)
                    return (df_out, 'tumor samples only (transformed)')
                elif force_raw == True:
                    df_out = df_targets_tumors[df_targets_tumors.TumorStageStr.isin(stage_choices)].copy()
                    #self.UpdateStagesList(df_out)
                    return (df_out, 'tumor samples only (untransformed)')
            elif dataTransformed == False:
                df_out = df_targets_tumors[df_targets_tumors.TumorStageStr.isin(stage_choices)].copy()
                #self.UpdateStagesList(df_out)
                return (df_out, 'tumor samples only (untransformed)')

    def CalculateCorrelation(self):
        TARGET_GENE = str(self.loaded_targets_list.currentItem().text())
        df_temp, samp_choice = get_sample_choice_df(self, force_raw=True)
        self.correlation_output_box.clear()
        data_summary_header = '-'*25+ 'Metadata for Analyzed Samples' + '-'*25
        self.correlation_output_box.appendPlainText(data_summary_header)

        data_summary_header = 'Tissue Type:\t\tNumber of Samples:'
        self.correlation_output_box.appendPlainText(data_summary_header)
        sample_types_dict = dict(df_temp.SampleType.value_counts())
        cum_sum = 0
        for key in sample_types_dict:
            if len(key.split(' ')) >= 2:
                msg='     -%s\t%d' % (key, sample_types_dict[key])
            else:
                msg='     -%s\t\t%d' % (key, sample_types_dict[key])
            cum_sum = cum_sum + int(sample_types_dict[key])
            self.correlation_output_box.appendPlainText(msg)

        data_summary_header = 'Tumor Stage:\t\tNumber of Samples'
        self.correlation_output_box.appendPlainText(data_summary_header)
        tumor_stages_dict = dict(df_temp.TumorStageStr.value_counts())
        cum_sum = 0
        for key in tumor_stages_dict:
            msg='     -Stage %s\t\t%d' % (key, tumor_stages_dict[key])
            cum_sum = cum_sum + int(tumor_stages_dict[key])
            self.correlation_output_box.appendPlainText(msg)
        msg='\tTotal:\t%d' % cum_sum
        self.correlation_output_box.appendPlainText(msg)
        msg = '-'*100
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
        global haveScatter, haveHeatmap, haveHistogram
        if ((self.loaded_targets_list.count() == 0) | (self.loaded_targets_list.currentItem() == None)):
            self.show_msg_box('Select a gene to plot', 'A gene must be selected from the list of loaded genes in order to generate histogram')
        else:
            TARGET_GENE = str(self.loaded_targets_list.currentItem().text())
            fig1 = Figure()
            ax1f1 = fig1.add_subplot(111)
            df_temp, samp_choice = get_sample_choice_df(self)
            print 'Loaded %s' % samp_choice
            ax1f1.hist(df_temp[TARGET_GENE])
            self.addmpl(fig1)
            haveScatter, haveHeatmap, haveHistogram = False, False, True
            print TARGET_GENE

    def GenerateScatterPlot(self):
        global haveScatter, haveHeatmap, haveHistogram
        if self.loaded_targets_list.count() == 0:
            self.show_msg_box('Select two genes to plot', 'Two genes must be selected from the list of loaded genes in order to generate scatterplot')
        else:
            selected_genes = self.loaded_targets_list.selectedItems()
            if len(selected_genes) == 2:
                gene_x = str(selected_genes[0].text())
                gene_y = str(selected_genes[1].text())
                df_temp, samp_choice = get_sample_choice_df(self)
                #if self.heatmap_windowed_checkbox.isChecked():
                #    sns.jointplot(x=gene_x, y=gene_y, data=df_temp, kind='reg', size=10, space=0)
                #    plt.show()
                #else:
                if self.save_figure_checkbox.isChecked():
                    fileName = str(QtGui.QFileDialog.getSaveFileName(self, 'Save as...', '/Users/TS_MBP/Documents/GitHub/TCGA-Tools-GUI/', selectedFilter='*.svg'))
                    file_format = fileName.split('.')[-1]
                    fig1 = Figure()
                    ax1f1 = fig1.add_subplot(111)
                    ax1f1.scatter(x=gene_x, y=gene_y, data=df_temp)
                    ax1f1.set(xlabel=gene_x, ylabel=gene_y)
                    self.addmpl(fig1)
                    fig1.savefig(fileName, format=file_format)
                    haveScatter, haveHeatmap, haveHistogram = True, False, True
                else:
                    fig1 = Figure()
                    ax1f1 = fig1.add_subplot(111)
                    ax1f1.scatter(x=gene_x, y=gene_y, data=df_temp)
                    ax1f1.set(xlabel=gene_x, ylabel=gene_y)
                    self.addmpl(fig1)
                    haveScatter, haveHeatmap, haveHistogram = True, False, True
            else:
                self.show_msg_box('Select two genes to plot', 'At least two genes must be selected from the list of loaded genes in order to generate scatterplot')

    def GenerateHeatMap(self):
        global haveScatter, haveHeatmap, haveHistogram
        if self.loaded_targets_list.count() == 0:
            self.show_msg_box('Select a target gene!', 'Select a gene from the list of loaded genes to serve as your gene of interest.  Target gene will appear first on the heatmap and samples (rows) will be sorted by target gene values.')
            return
        else:
            TARGET_GENE = str(self.loaded_targets_list.currentItem().text())
        df_temp, samp_choice = get_sample_choice_df(self)
        target_col = df_temp[TARGET_GENE]
        df_temp = df_temp.drop(labels=[TARGET_GENE], axis=1)
        df_temp.insert(0, TARGET_GENE, target_col)
        df_temp = df_temp.sort_values(TARGET_GENE, ascending=False)

        if self.heatmap_windowed_checkbox.isChecked():
            fig2 = Figure()
            ax1f2 = fig1.add_subplot(111)
            for item in (sns.heatmap(df_temp.select_dtypes(exclude=['object']),yticklabels=False, xticklabels=True, ax=ax1f2)).get_xticklabels():
                item.set_rotation(90)
            fig2.show()
        else:
            if self.save_figure_checkbox.isChecked():
                fileName = str(QtGui.QFileDialog.getSaveFileName(self, 'Save as...', '/Users/TS_MBP/Documents/GitHub/TCGA-Tools-GUI/', selectedFilter='*.svg'))
                file_format = fileName.split('.')[-1]
                fig1 = Figure()
                ax1f1 = fig1.add_subplot(111)
                haveScatter, haveHeatmap, haveHistogram = False, True, True
                for item in (sns.heatmap(df_temp.select_dtypes(exclude=['object']),yticklabels=False, xticklabels=True, ax=ax1f1)).get_xticklabels():
                    item.set_rotation(90)
                fig1.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.12)
                self.addmpl(fig1)
                fig1.savefig(fileName, format=file_format)
                haveScatter, haveHeatmap, haveHistogram = False, True, True
            else:
                fig1 = Figure()
                ax1f1 = fig1.add_subplot(111)
                for item in (sns.heatmap(df_temp.select_dtypes(exclude=['object']),yticklabels=False, xticklabels=True, ax=ax1f1)).get_xticklabels():
                    item.set_rotation(90)
                fig1.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.12)
                self.addmpl(fig1)
                haveScatter, haveHeatmap, haveHistogram = False, True, True


    def ExportToExcel(self):
        fileName = str(QtGui.QFileDialog.getSaveFileName(self, 'Save as...', '/Users/TS_MBP/Documents/GitHub/TCGA-Tools-GUI/', selectedFilter='*.xlsx'))
        df_temp, samp_choice = get_sample_choice_df(self)
        baseFileName = fileName.split('.')[0]
        excelFileName = '%s.xlsx' % baseFileName
        df_temp.to_excel(excelFileName, index=False)
        self.log_display_box.appendPlainText(("Saved %s data to Excel table as %s" % (samp_choice, excelFileName)))

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    w = QWidget()
    window = MyApp()
    window.show()
    sys.exit(app.exec_())

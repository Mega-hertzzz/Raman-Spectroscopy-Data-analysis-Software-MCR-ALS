# -*- coding: utf-8 -*-
# Form implementation generated from reading ui file 'Radans2.ui'
# Creator : Mohammad Hosein Zakaryapour (MHz: MegaHertz)

# PARS: a Program for Analysis of Raman Spectra


import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
import pandas as pd
import rampy as rp
import scipy.signal as ss
import os 
import fnmatch
import random 
import statistics as st
import sys
import time
import h5py
import math
import despike
import scipy.io
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sklearn.cluster import KMeans
from numpy.linalg import inv
from sklearn import decomposition
from sympy import pretty_print as latex
from pymcr.mcr import McrAR
from pymcr.regressors import OLS, NNLS
from pymcr.constraints import ConstraintNonneg, ConstraintNorm
from pymcr.metrics import mse
from pymcr.constraints import Constraint
from scipy import signal
from lmfit.models import GaussianModel, ConstantModel
from sklearn.linear_model import Ridge, Lasso, LinearRegression
from sklearn.preprocessing import minmax_scale
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog, QMessageBox , QPushButton, QVBoxLayout, QDialog
from PyQt5.QtGui import QIcon
from PyQt5.Qt import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar



class Ui_MainWindow(QWidget): #(object)

    
    def setupUi(self, MainWindow, parent=None):
        super(Ui_MainWindow, self).__init__(parent)

        self.check_open = 0
        self.check_save_preprocessing = 0
        self.check_pca = 0
        self.check_colorbar = 0
        self.check_clsref1 = 0
        self.check_clsref2 = 0
        self.check_clssam = 0
        self.check_kmeans = 0
        self.check_kmeans_clear = 0
        self.check_run_preprocessing = 0
        self.norm_snv_check = 0
        self.norm_area_check = 0
        self.check_save_kmeans = 0
        self.check_mcr_clear = 0
        self.check_kmca = 0
        self.check_nmf = 0
        self.check_nmf_clear = 0
        self.check_save_mcr = 0
        self.check_save_nmf = 0
        self.norm_area_check_kmeans = 0
        self.laser_wavelength = 780
        self.note_unit_conversion = " " 
        self.note_smoothing = " "
        self.note_baseline = " "
        self.note_normalization = " "

      

        scriptDir = os.path.dirname(os.path.realpath(__file__))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(scriptDir + os.path.sep + "ice-cream.png"), QtGui.QIcon.Selected, QtGui.QIcon.On)
        MainWindow.setWindowIcon(icon)

        MainWindow.setObjectName("MainWindow")
        # MainWindow.resize(1047, 650)
        MainWindow.setGeometry(60,60,1047, 650)

        os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"
        # qapp = QApplication(sys.argv)
        # qapp.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)

        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        MainWindow.setFont(font)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(0, 0, 1045, 650))
        
        ### Deleted the usless stuff


        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.tabWidget.setFont(font)
        self.tabWidget.setStatusTip("")
        self.tabWidget.setStyleSheet("background-color: rgb(240, 240, 240);\n"
"selection-background-color: rgb(255, 255, 255);")
        self.tabWidget.setObjectName("tabWidget")
        self.tab1_preprocessing = QtWidgets.QWidget()
        self.tab1_preprocessing.setObjectName("tab1_preprocessing")
        self.frame = QtWidgets.QFrame(self.tab1_preprocessing)
        self.frame.setGeometry(QtCore.QRect(-1, 0, 1041, 101))  

        ### Deleted the usless stuff

        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.frame.setFont(font)
        self.frame.setFrameShape(QtWidgets.QFrame.Box)
        self.frame.setFrameShadow(QtWidgets.QFrame.Plain)
        self.frame.setObjectName("frame")
        self.gridLayout = QtWidgets.QGridLayout(self.frame)
        self.gridLayout.setObjectName("gridLayout")
        self.smotthing_box = QtWidgets.QGroupBox(self.frame)

        ### Deleted the usless stuff ###

        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(11)
        self.smotthing_box.setFont(font)
        self.smotthing_box.setObjectName("smotthing_box")
        self.lineEdit_4_OrderOfPoly = QtWidgets.QLineEdit(self.smotthing_box)
        self.lineEdit_4_OrderOfPoly.setGeometry(QtCore.QRect(130, 50, 21, 20))

        ### Deleted the usless stuff  ...
        palette = QtGui.QPalette()
        
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.HighlightedText, brush)
        
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.HighlightedText, brush)
        
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.HighlightedText, brush)
        self.lineEdit_4_OrderOfPoly.setPalette(palette)
        ### ...

        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.lineEdit_4_OrderOfPoly.setFont(font)
        self.lineEdit_4_OrderOfPoly.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.lineEdit_4_OrderOfPoly.setAlignment(QtCore.Qt.AlignCenter)
        self.lineEdit_4_OrderOfPoly.setObjectName("lineEdit_4_OrderOfPoly")
        self.label_9_order_of_poly = QtWidgets.QLabel(self.smotthing_box)
        self.label_9_order_of_poly.setGeometry(QtCore.QRect(10, 50, 119, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        font.setBold(False)
        font.setWeight(50)
        self.label_9_order_of_poly.setFont(font)
        self.label_9_order_of_poly.setObjectName("label_9_order_of_poly")
        self.lineEdit_3_WindowLength = QtWidgets.QLineEdit(self.smotthing_box)
        self.lineEdit_3_WindowLength.setGeometry(QtCore.QRect(130, 22, 21, 20))

        ### Deleted the usless stuff ___
        palette = QtGui.QPalette()
        
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.HighlightedText, brush)
        
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.HighlightedText, brush)
        
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.HighlightedText, brush)
        self.lineEdit_3_WindowLength.setPalette(palette)
        ### ___

        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.lineEdit_3_WindowLength.setFont(font)
        self.lineEdit_3_WindowLength.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.lineEdit_3_WindowLength.setAlignment(QtCore.Qt.AlignCenter)
        self.lineEdit_3_WindowLength.setObjectName("lineEdit_3_WindowLength")
        self.label_10_window_length = QtWidgets.QLabel(self.smotthing_box)
        self.label_10_window_length.setGeometry(QtCore.QRect(10, 24, 111, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        font.setBold(False)
        font.setWeight(50)
        self.label_10_window_length.setFont(font)
        self.label_10_window_length.setObjectName("label_0_window_length")

        self.label_20_window_length = QtWidgets.QLabel(self.smotthing_box)
        self.label_20_window_length.setGeometry(QtCore.QRect(190, 50, 100, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        font.setBold(False)
        font.setWeight(50)
        self.label_20_window_length.setFont(font)
        self.label_20_window_length.setObjectName("label_20_window_length")
        self.label_20_window_length.setEnabled(False)

        self.lineEdit_9_wl_despiking = QtWidgets.QLineEdit(self.smotthing_box)
        self.lineEdit_9_wl_despiking.setGeometry(QtCore.QRect(259, 50, 34, 20))
        self.lineEdit_9_wl_despiking.setFont(font)
        self.lineEdit_9_wl_despiking.setAutoFillBackground(False)
        self.lineEdit_9_wl_despiking.setAlignment(QtCore.Qt.AlignCenter)
        self.lineEdit_9_wl_despiking.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.lineEdit_9_wl_despiking.setText("")
        self.lineEdit_9_wl_despiking.setObjectName("lineEdit_9_wl_despiking")
        self.lineEdit_9_wl_despiking.setPalette(palette)
        self.lineEdit_9_wl_despiking.setEnabled(False)

        self.checkbox_despiking = QtWidgets.QCheckBox(self.smotthing_box)
        self.checkbox_despiking.setGeometry(QtCore.QRect(190, 28, 110, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        font.setStyleStrategy(QtGui.QFont.PreferDefault)
        self.checkbox_despiking.setFont(font)
        self.checkbox_despiking.setObjectName("checkbox_despiking")

        self.gridLayout.addWidget(self.smotthing_box, 0, 2, 5, 1)
        self.Button_save = QtWidgets.QPushButton(self.frame)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_save.setFont(font)
        self.Button_save.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_save.setObjectName("Button_save")
        self.gridLayout.addWidget(self.Button_save, 2, 0, 1, 1)
        self.spectral_range_box = QtWidgets.QGroupBox(self.frame)
        
        ### Deleted the usless stuff ###

        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(11)
        self.spectral_range_box.setFont(font)
        self.spectral_range_box.setObjectName("spectral_range_box")
        self.label_5_cm = QtWidgets.QLabel(self.spectral_range_box)
        self.label_5_cm.setGeometry(QtCore.QRect(85, 49, 41, 21))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(9)
        self.label_5_cm.setFont(font)
        self.label_5_cm.setObjectName("label_5_cm")
        self.lineEdit_1_range = QtWidgets.QLineEdit(self.spectral_range_box)
        self.lineEdit_1_range.setGeometry(QtCore.QRect(45, 50, 36, 20))

        ### Deleted the usless stuff ___
        palette = QtGui.QPalette()
        
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.HighlightedText, brush)
        
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.HighlightedText, brush)
        
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.HighlightedText, brush)
        self.lineEdit_1_range.setPalette(palette)
        ### ___

        font = QtGui.QFont()
        font.setPointSize(10)
        self.lineEdit_1_range.setFont(font)
        self.lineEdit_1_range.setAutoFillBackground(False)
        self.lineEdit_1_range.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.lineEdit_1_range.setText("")
        self.lineEdit_1_range.setObjectName("lineEdit_1_range")
        self.label_6_cm = QtWidgets.QLabel(self.spectral_range_box)
        self.label_6_cm.setGeometry(QtCore.QRect(180, 50, 41, 21))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(9)
        self.label_6_cm.setFont(font)
        self.label_6_cm.setObjectName("label_6_cm")


        self.label_12_nm = QtWidgets.QLabel(self.spectral_range_box)
        self.label_12_nm.setGeometry(QtCore.QRect(280, 25, 20, 21))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(9)
        self.label_12_nm.setFont(font)
        self.label_12_nm.setObjectName("label_12_nm")
        self.label_12_nm.setEnabled(False)

        self.lineEdit_5_laser = QtWidgets.QLineEdit(self.spectral_range_box)
        self.lineEdit_5_laser.setGeometry(QtCore.QRect(239, 24, 36, 20))
        self.lineEdit_5_laser.setPalette(palette)
        self.lineEdit_5_laser.setFont(font)
        self.lineEdit_5_laser.setAutoFillBackground(False)
        self.lineEdit_5_laser.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.lineEdit_5_laser.setText("")
        self.lineEdit_5_laser.setObjectName("lineEdit_5_laser")
        self.lineEdit_5_laser.setEnabled(False)

        self.lineEdit_2_range = QtWidgets.QLineEdit(self.spectral_range_box)
        self.lineEdit_2_range.setGeometry(QtCore.QRect(138, 50, 36, 20))

        ### Deleted the usless stuff ...
        palette = QtGui.QPalette()
        
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.HighlightedText, brush)
        
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.HighlightedText, brush)
        
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.HighlightedText, brush)
        self.lineEdit_2_range.setPalette(palette)
        ### ...

        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        self.lineEdit_2_range.setFont(font)
        self.lineEdit_2_range.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.lineEdit_2_range.setText("")
        self.lineEdit_2_range.setObjectName("lineEdit_2_range")
        self.label_7_to = QtWidgets.QLabel(self.spectral_range_box)
        self.label_7_to.setGeometry(QtCore.QRect(121, 52, 16, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        self.label_7_to.setFont(font)
        self.label_7_to.setObjectName("label_7_to")
        self.checkbox_change_unit = QtWidgets.QCheckBox(self.spectral_range_box)
        self.checkbox_change_unit.setGeometry(QtCore.QRect(10, 25, 211, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        font.setStyleStrategy(QtGui.QFont.PreferDefault)
        self.checkbox_change_unit.setFont(font)
        self.checkbox_change_unit.setObjectName("checkbox_change_unit")
        self.label_8_from = QtWidgets.QLabel(self.spectral_range_box)
        self.label_8_from.setGeometry(QtCore.QRect(10, 51, 31, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        self.label_8_from.setFont(font)
        self.label_8_from.setObjectName("label_8_from")

        self.label_11_laser = QtWidgets.QLabel(self.spectral_range_box)
        self.label_11_laser.setGeometry(QtCore.QRect(200, 25, 35, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        self.label_11_laser.setFont(font)
        self.label_11_laser.setObjectName("label_11_laser")
        self.label_11_laser.setEnabled(False)


        self.gridLayout.addWidget(self.spectral_range_box, 0, 1, 5, 1)
        self.Button_open = QtWidgets.QPushButton(self.frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Button_open.sizePolicy().hasHeightForWidth())
        self.Button_open.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_open.setFont(font)
        self.Button_open.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_open.setStyleSheet("border-top-color: rgb(255, 1, 81);\n"
"border-top-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 rgba(0, 0, 0, 255), stop:1 rgba(255, 255, 255, 255));")
        self.Button_open.setObjectName("Button_open")
        self.gridLayout.addWidget(self.Button_open, 0, 0, 2, 1)
        self.Button_run1_preprocessing = QtWidgets.QPushButton(self.frame)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_run1_preprocessing.setFont(font)
        self.Button_run1_preprocessing.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_run1_preprocessing.setObjectName("Button_run1_preprocessing")
        self.gridLayout.addWidget(self.Button_run1_preprocessing, 3, 0, 2, 1)

        self.Others_preprocessing_box = QtWidgets.QGroupBox(self.frame)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(11)
        self.Others_preprocessing_box.setFont(font)
        self.Others_preprocessing_box.setTitle("")
        self.Others_preprocessing_box.setFlat(False)
        self.Others_preprocessing_box.setObjectName("Others_preprocessing_box")
        self.frame_baseline = QtWidgets.QFrame(self.Others_preprocessing_box)
        self.frame_baseline.setGeometry(QtCore.QRect(6, 1, 301, 41))
        self.frame_baseline.setToolTipDuration(0)
        self.frame_baseline.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_baseline.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_baseline.setObjectName("frame_baseline")
        self.label_baseline = QtWidgets.QLabel(self.frame_baseline)
        self.label_baseline.setGeometry(QtCore.QRect(3, -1, 201, 20))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(11)
        self.label_baseline.setFont(font)
        self.label_baseline.setToolTipDuration(0)
        self.label_baseline.setObjectName("label_baseline")
        self.radioButton_baseline_arpls = QtWidgets.QRadioButton(self.frame_baseline)
        self.radioButton_baseline_arpls.setGeometry(QtCore.QRect(128, 17, 61, 17))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        self.radioButton_baseline_arpls.setFont(font)
        self.radioButton_baseline_arpls.setChecked(False)
        self.radioButton_baseline_arpls.setObjectName("radioButton_baseline_arpls")
        self.radioButton_baseline_drpls = QtWidgets.QRadioButton(self.frame_baseline)
        self.radioButton_baseline_drpls.setGeometry(QtCore.QRect(200, 17, 61, 17))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        self.radioButton_baseline_drpls.setFont(font)
        self.radioButton_baseline_drpls.setChecked(False)
        self.radioButton_baseline_drpls.setObjectName("radioButton_baseline_drpls")

        self.radioButton_baseline_none = QtWidgets.QRadioButton(self.frame_baseline)
        self.radioButton_baseline_none.setGeometry(QtCore.QRect(56, 17, 61, 17))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        self.radioButton_baseline_none.setFont(font)
        self.radioButton_baseline_none.setChecked(True)
        self.radioButton_baseline_none.setObjectName("radioButton_baseline_drpls")

        self.frame_normalization = QtWidgets.QFrame(self.Others_preprocessing_box)
        self.frame_normalization.setGeometry(QtCore.QRect(5, 45, 301, 31))
        self.frame_normalization.setToolTipDuration(0)
        self.frame_normalization.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_normalization.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_normalization.setObjectName("frame_normalization")
        self.label_normalization = QtWidgets.QLabel(self.frame_normalization)
        self.label_normalization.setGeometry(QtCore.QRect(5, 1, 171, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(11)
        self.label_normalization.setFont(font)
        self.label_normalization.setToolTipDuration(0)
        self.label_normalization.setObjectName("label_normalization")
        self.radioButton_normalization_snv = QtWidgets.QRadioButton(self.frame_normalization)
        self.radioButton_normalization_snv.setGeometry(QtCore.QRect(161, 15, 51, 20))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        self.radioButton_normalization_snv.setFont(font)
        self.radioButton_normalization_snv.setChecked(False)
        self.radioButton_normalization_snv.setObjectName("radioButton_normalization_snv")
        self.radioButton_normalization_area = QtWidgets.QRadioButton(self.frame_normalization)
        self.radioButton_normalization_area.setGeometry(QtCore.QRect(89, 15, 51, 20))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        self.radioButton_normalization_area.setFont(font)
        self.radioButton_normalization_area.setChecked(True)
        self.radioButton_normalization_area.setObjectName("radioButton_normalization_area")
        self.gridLayout.addWidget(self.Others_preprocessing_box, 0, 3, 5, 1)
        self.frame_3 = QtWidgets.QFrame(self.tab1_preprocessing)
        self.frame_3.setGeometry(QtCore.QRect(-8, 99, 1050, 511))
        font = QtGui.QFont()
        font.setFamily("Roboto")

        self.layout_3 = QtWidgets.QGridLayout(self.frame_3) ######### added for the plot

        self.frame_3.setFont(font)
        self.frame_3.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.frame_3.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_3.setObjectName("frame_3")
        self.frame_3.raise_()
        self.frame.raise_()
        self.tabWidget.addTab(self.tab1_preprocessing, "")
        self.tab_2_pca = QtWidgets.QWidget()
        self.tab_2_pca.setObjectName("tab_2_pca")
        self.tabWidget_2 = QtWidgets.QTabWidget(self.tab_2_pca)
        self.tabWidget_2.setGeometry(QtCore.QRect(-4, 29, 1051, 601))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.tabWidget_2.setFont(font)
        self.tabWidget_2.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.tabWidget_2.setObjectName("tabWidget_2")
        self.tab_3_scatter = QtWidgets.QWidget()
        self.tab_3_scatter.setObjectName("tab_3_scatter")
        self.tabWidget_2.addTab(self.tab_3_scatter, "")
        self.tab_4_loading = QtWidgets.QWidget()
        self.tab_4_loading.setObjectName("tab_4_loading")
        self.tabWidget_2.addTab(self.tab_4_loading, "")
        self.tab_6_pca_image = QtWidgets.QWidget()
        self.tab_6_pca_image.setObjectName("tab_6_pca_image")
        self.tabWidget_2.addTab(self.tab_6_pca_image, "")
        self.tab_5_scree = QtWidgets.QWidget()
        self.tab_5_scree.setObjectName("tab_5_scree")
        self.tabWidget_2.addTab(self.tab_5_scree, "")
        

        

        ######## Made frames for pca plots ########

        self.frame_4 = QtWidgets.QFrame(self.tab_3_scatter)
        self.frame_4.setGeometry(QtCore.QRect(250, 0, 550, 550))
        self.frame_4.setStyleSheet("background-color: rgb(255, 255, 255);")

        self.frame_5 = QtWidgets.QFrame(self.tab_4_loading)
        self.frame_5.setGeometry(QtCore.QRect(-7, 0, 1041, 550))
        self.frame_5.setStyleSheet("background-color: rgb(255, 255, 255);")

        self.frame_6 = QtWidgets.QFrame(self.tab_5_scree)
        self.frame_6.setGeometry(QtCore.QRect(250, 0, 550, 550))
        self.frame_6.setStyleSheet("background-color: rgb(255, 255, 255);")

        self.frame_7 = QtWidgets.QFrame(self.tab_6_pca_image)
        self.frame_7.setGeometry(QtCore.QRect(-7, 0, 1041, 550))
        self.frame_7.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.frame_7.setObjectName("frame_7")
        self.frame_7.raise_()


        
        self.layout_4 = QtWidgets.QGridLayout(self.frame_4) ######### added for the plot
        self.layout_5 = QtWidgets.QGridLayout(self.frame_5)
        self.layout_6 = QtWidgets.QGridLayout(self.frame_6)
        self.layout_7 = QtWidgets.QGridLayout(self.frame_7) #########


        # The secret button for changing the colormap of the images inside pca score image tab
        self.radioButton_pcacolormap = QtWidgets.QRadioButton(self.frame_7)
        self.radioButton_pcacolormap.setGeometry(QtCore.QRect(1018, 536, 51, 20))
        self.radioButton_pcacolormap.setChecked(False)
        self.radioButton_pcacolormap.setObjectName("radioButton_pcacolormap")

        # Kmeans for PCA checkbox
        self.checkbox_kmean_pca = QtWidgets.QCheckBox(self.tab_3_scatter)
        self.checkbox_kmean_pca.setGeometry(QtCore.QRect(271, 50, 300, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        font.setStyleStrategy(QtGui.QFont.PreferDefault)
        self.checkbox_kmean_pca.setFont(font)
        self.checkbox_kmean_pca.setObjectName("checkbox_kmean_pca")

        # Label for Kmeans for PCA checkbox
        self.label_13_pca_clusters = QtWidgets.QLabel(self.tab_3_scatter)
        self.label_13_pca_clusters.setGeometry(QtCore.QRect(271, 74, 200, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        self.label_13_pca_clusters.setFont(font)
        self.label_13_pca_clusters.setObjectName("label_13_pca_clusters")
        self.label_13_pca_clusters.setEnabled(False)

        # Lineedit for Kmeans for PCA checkbox
        self.lineEdit_6_Kmeans_pca = QtWidgets.QLineEdit(self.tab_3_scatter)
        self.lineEdit_6_Kmeans_pca.setGeometry(QtCore.QRect(390, 73, 30, 20))
        self.lineEdit_6_Kmeans_pca.setPalette(palette)
        self.lineEdit_6_Kmeans_pca.setFont(font)
        self.lineEdit_6_Kmeans_pca.setAutoFillBackground(False)
        self.lineEdit_6_Kmeans_pca.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.lineEdit_6_Kmeans_pca.setText("")
        self.lineEdit_6_Kmeans_pca.setObjectName("lineEdit_6_Kmeans_pca")
        self.lineEdit_6_Kmeans_pca.setEnabled(False)

        self.frame_2 = QtWidgets.QFrame(self.tab_2_pca)
        self.frame_2.setGeometry(QtCore.QRect(0, 0, 1041, 51))
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")

        self.Button_run2_pca = QtWidgets.QPushButton(self.frame_2)
        self.Button_run2_pca.setGeometry(QtCore.QRect(480, 6, 75, 38))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Button_run2_pca.sizePolicy().hasHeightForWidth())
        self.Button_run2_pca.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_run2_pca.setFont(font)
        self.Button_run2_pca.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_run2_pca.setObjectName("Button_run2_pca")
        self.label = QtWidgets.QLabel(self.frame_2)
        self.label.setGeometry(QtCore.QRect(730, 4, 351, 41))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.frame_2.raise_()
        self.tabWidget_2.raise_()
        self.tabWidget.addTab(self.tab_2_pca, "")
        self.tab_3_cls = QtWidgets.QWidget()
        self.tab_3_cls.setObjectName("tab_3_cls")

        self.frame_8 = QtWidgets.QFrame(self.tab_3_cls)
        self.frame_8.setGeometry(QtCore.QRect(-8, 99, 1050, 511))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.frame_8.setFont(font)
        self.frame_8.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.frame_8.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_8.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_8.setObjectName("frame_8")
        self.frame_5_cls = QtWidgets.QFrame(self.tab_3_cls)
        self.frame_5_cls.setGeometry(QtCore.QRect(-1, 0, 1041, 101))

        self.layout_8 = QtWidgets.QGridLayout(self.frame_8) ###### added for the cls plot
              
        self.tab_4_kmeans = QtWidgets.QWidget()
        self.tab_4_kmeans.setObjectName("tab_4_kmeans")

        self.tab_5_mcr = QtWidgets.QWidget()
        self.tab_5_mcr.setObjectName("tab_5_mcr")

        self.tab_6_nmf = QtWidgets.QWidget()
        self.tab_6_nmf.setObjectName("tab_6_nmf")

        ### Deleted the usless stuff ###

        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.frame_5_cls.setFont(font)
        self.frame_5_cls.setFrameShape(QtWidgets.QFrame.Box)
        self.frame_5_cls.setFrameShadow(QtWidgets.QFrame.Plain)
        self.frame_5_cls.setObjectName("frame_5_cls")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.frame_5_cls)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.Button_reference_1 = QtWidgets.QPushButton(self.frame_5_cls)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Button_reference_1.sizePolicy().hasHeightForWidth())
        self.Button_reference_1.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_reference_1.setFont(font)
        self.Button_reference_1.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_reference_1.setStyleSheet("border-top-color: rgb(255, 1, 81);\n"
"border-top-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 rgba(0, 0, 0, 255), stop:1 rgba(255, 255, 255, 255));")
        self.Button_reference_1.setObjectName("Button_reference_1")
        self.gridLayout_2.addWidget(self.Button_reference_1, 2, 0, 1, 1)
        self.frame_6_cls = QtWidgets.QFrame(self.frame_5_cls)
        self.frame_6_cls.setFrameShape(QtWidgets.QFrame.Box)
        self.frame_6_cls.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_6_cls.setObjectName("frame_6_cls")
        self.lcdNumber_component1 = QtWidgets.QLCDNumber(self.frame_6_cls)
        self.lcdNumber_component1.setGeometry(QtCore.QRect(360, 9, 141, 31))
        self.lcdNumber_component1.setDigitCount(10)
        self.lcdNumber_component1.setObjectName("lcdNumber_component1")
        self.lcdNumber_component1.setMode(QtWidgets.QLCDNumber.Dec)
        self.lcdNumber_component1.setSegmentStyle(QtWidgets.QLCDNumber.Flat)
        self.lcdNumber_component2 = QtWidgets.QLCDNumber(self.frame_6_cls)
        self.lcdNumber_component2.setGeometry(QtCore.QRect(360, 41, 141, 31))
        self.lcdNumber_component2.setSmallDecimalPoint(False)
        self.lcdNumber_component2.setDigitCount(10)
        self.lcdNumber_component2.setMode(QtWidgets.QLCDNumber.Dec)
        self.lcdNumber_component2.setSegmentStyle(QtWidgets.QLCDNumber.Flat)
        self.lcdNumber_component2.setObjectName("lcdNumber_component2")
        self.label_2_ratio1 = QtWidgets.QLabel(self.frame_6_cls)
        self.label_2_ratio1.setGeometry(QtCore.QRect(180, 16, 180, 20))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2_ratio1.sizePolicy().hasHeightForWidth())
        self.label_2_ratio1.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(11)
        self.label_2_ratio1.setFont(font)
        self.label_2_ratio1.setObjectName("label_2_ratio1")
        self.label_3_ratio2 = QtWidgets.QLabel(self.frame_6_cls)
        self.label_3_ratio2.setGeometry(QtCore.QRect(180, 47, 180, 20))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(11)
        self.label_3_ratio2.setFont(font)
        self.label_3_ratio2.setObjectName("label_3_ratio2")

        

        ###### Kmeans Clustering #########################
        #
        #
        self.frame_9 = QtWidgets.QFrame(self.tab_4_kmeans)
        self.frame_9.setGeometry(QtCore.QRect(0, 0, 1041, 51))
        self.frame_9.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_9.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_9.setObjectName("frame_9")
        self.frame_9.raise_()

        self.Button_run_kmeans = QtWidgets.QPushButton(self.frame_9)
        self.Button_run_kmeans.setGeometry(QtCore.QRect(520, 6, 75, 38))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_run_kmeans.setFont(font)
        self.Button_run_kmeans.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_run_kmeans.setObjectName("Button_run_kmeans")

        self.Button_kmeans_save = QtWidgets.QPushButton(self.frame_9)
        self.Button_kmeans_save.setGeometry(QtCore.QRect(440, 6, 75, 38))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_kmeans_save.setFont(font)
        self.Button_kmeans_save.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_kmeans_save.setObjectName("Button_kmeans_save")

        self.tabWidget_kmeans = QtWidgets.QTabWidget(self.tab_4_kmeans)
        self.tabWidget_kmeans.setGeometry(QtCore.QRect(-4, 29, 1051, 601))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.tabWidget_kmeans.setFont(font)
        self.tabWidget_kmeans.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.tabWidget_kmeans.setObjectName("tabWidget_kmeans")
        self.tab_8_kmeans_mean = QtWidgets.QWidget()
        self.tab_8_kmeans_mean.setObjectName("tab_8_kmeans_mean")
        self.tabWidget_kmeans.addTab(self.tab_8_kmeans_mean, "")
        self.tab_7_kmeans_image = QtWidgets.QWidget()
        self.tab_7_kmeans_image.setObjectName("tab_7_kmeans_image")
        self.tabWidget_kmeans.addTab(self.tab_7_kmeans_image, "")
        self.tab_9_elbow = QtWidgets.QWidget()
        self.tab_9_elbow.setObjectName("tab_9_elbow")
        self.tabWidget_kmeans.addTab(self.tab_9_elbow, "")


        self.label_14_kmeans_num = QtWidgets.QLabel(self.frame_9)
        self.label_14_kmeans_num.setGeometry(QtCore.QRect(125, 6, 150, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        self.label_14_kmeans_num.setFont(font)
        self.label_14_kmeans_num.setObjectName("label_14_kmeans_num")

        self.lineEdit_7_kmeans_clusters = QtWidgets.QLineEdit(self.frame_9)
        self.lineEdit_7_kmeans_clusters.setGeometry(QtCore.QRect(245, 5, 30, 20))
        self.lineEdit_7_kmeans_clusters.setPalette(palette)
        self.lineEdit_7_kmeans_clusters.setFont(font)
        self.lineEdit_7_kmeans_clusters.setAutoFillBackground(False)
        self.lineEdit_7_kmeans_clusters.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.lineEdit_7_kmeans_clusters.setText("")
        self.lineEdit_7_kmeans_clusters.setObjectName("lineEdit_7_kmeans_clusters")
        
        self.frame_10 = QtWidgets.QFrame(self.tab_7_kmeans_image)
        self.frame_10.setGeometry(QtCore.QRect(-7, 0, 1041, 550))
        self.frame_10.setStyleSheet("background-color: rgb(255, 255, 255);")

        self.frame_11 = QtWidgets.QFrame(self.tab_8_kmeans_mean)
        self.frame_11.setGeometry(QtCore.QRect(-7, 0, 1041, 560))
        self.frame_11.setStyleSheet("background-color: rgb(255, 255, 255);")

        self.frame_12 = QtWidgets.QFrame(self.tab_9_elbow)
        self.frame_12.setGeometry(QtCore.QRect(-7, 0, 1041, 550))
        self.frame_12.setStyleSheet("background-color: rgb(255, 255, 255);")

        self.layout_9 = QtWidgets.QGridLayout(self.frame_10)
        self.layout_10 = QtWidgets.QGridLayout(self.frame_11)
        self.layout_11 = QtWidgets.QGridLayout(self.frame_12)

        self.label_15 = QtWidgets.QLabel(self.frame_9)
        self.label_15.setGeometry(QtCore.QRect(730, 4, 351, 41))
        sizePolicy.setHeightForWidth(self.label_15.sizePolicy().hasHeightForWidth())
        self.label_15.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.label_15.setFont(font)
        self.label_15.setObjectName("label")

        #
        #
        ###### End of Kmeans #############################################

        ###### MCR #######################################################
        #
        #
        self.frame_13 = QtWidgets.QFrame(self.tab_5_mcr)
        self.frame_13.setGeometry(QtCore.QRect(0, 0, 1041, 51))
        self.frame_13.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_13.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_13.setObjectName("frame_13")
        self.frame_13.raise_()

        self.Button_run_mcr = QtWidgets.QPushButton(self.frame_13)
        self.Button_run_mcr.setGeometry(QtCore.QRect(520, 6, 75, 38))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_run_mcr.setFont(font)
        self.Button_run_mcr.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_run_mcr.setObjectName("Button_run_mcr")

        self.Button_mcr_save = QtWidgets.QPushButton(self.frame_13)
        self.Button_mcr_save.setGeometry(QtCore.QRect(440, 6, 75, 38))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_mcr_save.setFont(font)
        self.Button_mcr_save.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_mcr_save.setObjectName("Button_mcr_save")

        self.tabWidget_mcr = QtWidgets.QTabWidget(self.tab_5_mcr)
        self.tabWidget_mcr.setGeometry(QtCore.QRect(-4, 29, 1051, 601))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.tabWidget_mcr.setFont(font)
        self.tabWidget_mcr.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.tabWidget_mcr.setObjectName("tabWidget_mcr")
        self.tab_9_mcr_spectra = QtWidgets.QWidget()
        self.tab_9_mcr_spectra.setObjectName("tab_9_mcr_spectra")
        self.tabWidget_mcr.addTab(self.tab_9_mcr_spectra, "")
        self.tab_10_mcr_image = QtWidgets.QWidget()
        self.tab_10_mcr_image.setObjectName("tab_10_mcr_image")
        self.tabWidget_mcr.addTab(self.tab_10_mcr_image, "")

        
        self.frame_14 = QtWidgets.QFrame(self.tab_9_mcr_spectra)
        self.frame_14.setGeometry(QtCore.QRect(-7, 0, 1041, 550))
        self.frame_14.setStyleSheet("background-color: rgb(255, 255, 255);")

        self.frame_15 = QtWidgets.QFrame(self.tab_10_mcr_image)
        self.frame_15.setGeometry(QtCore.QRect(-7, 0, 1055, 550))
        self.frame_15.setStyleSheet("background-color: rgb(255, 255, 255);")

        self.layout_12 = QtWidgets.QGridLayout(self.frame_14)
        self.layout_13 = QtWidgets.QGridLayout(self.frame_15)

        self.label_16 = QtWidgets.QLabel(self.frame_13)
        self.label_16.setGeometry(QtCore.QRect(730, 4, 351, 41))
        sizePolicy.setHeightForWidth(self.label_16.sizePolicy().hasHeightForWidth())
        self.label_16.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.label_16.setFont(font)
        self.label_16.setObjectName("label")

        self.frame_initial_spectra = QtWidgets.QFrame(self.frame_13)
        self.frame_initial_spectra.setGeometry(QtCore.QRect(1, 1, 400, 30))
        self.frame_initial_spectra.setToolTipDuration(0)
        self.frame_initial_spectra.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_initial_spectra.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_initial_spectra.setObjectName("frame_initial_spectra")

        self.radioButton_initial_nmf = QtWidgets.QRadioButton(self.frame_initial_spectra)
        self.radioButton_initial_nmf.setGeometry(QtCore.QRect(155, 5, 76, 17))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(9)
        self.radioButton_initial_nmf.setFont(font)
        self.radioButton_initial_nmf.setChecked(True)
        self.radioButton_initial_nmf.setObjectName("radioButton_initial_nmf")
        self.radioButton_initial_kmeans = QtWidgets.QRadioButton(self.frame_initial_spectra)
        self.radioButton_initial_kmeans.setGeometry(QtCore.QRect(240, 5, 61, 17))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(9)
        self.radioButton_initial_kmeans.setFont(font)
        self.radioButton_initial_kmeans.setChecked(False)
        self.radioButton_initial_kmeans.setObjectName("radioButton_initial_kmeans")

        self.label_19 = QtWidgets.QLabel(self.frame_initial_spectra)
        self.label_19.setGeometry(QtCore.QRect(57, 5, 90, 20))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        self.label_19.setFont(font)
        self.label_19.setToolTipDuration(0)
        self.label_19.setObjectName("label_19")


        #
        #
        ###### End of MCR  #############################################


        ###### NMF #######################################################
        #
        #
        self.frame_16 = QtWidgets.QFrame(self.tab_6_nmf)
        self.frame_16.setGeometry(QtCore.QRect(0, 0, 1041, 51))
        self.frame_16.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_16.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_16.setObjectName("frame_16")
        self.frame_16.raise_()

        self.Button_run_nmf = QtWidgets.QPushButton(self.frame_16)
        self.Button_run_nmf.setGeometry(QtCore.QRect(520, 6, 75, 38))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_run_nmf.setFont(font)
        self.Button_run_nmf.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_run_nmf.setObjectName("Button_run_mcr")

        self.Button_nmf_save = QtWidgets.QPushButton(self.frame_16)
        self.Button_nmf_save.setGeometry(QtCore.QRect(440, 6, 75, 38))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_nmf_save.setFont(font)
        self.Button_nmf_save.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_nmf_save.setObjectName("Button_mcr_save")

        
        self.frame_17 = QtWidgets.QFrame(self.tab_6_nmf)
        self.frame_17 .setGeometry(QtCore.QRect(-9, 51, 1050, 560))
        self.frame_17 .setStyleSheet("background-color: rgb(255, 255, 255);")

        self.label_17_nmf_num = QtWidgets.QLabel(self.frame_16)
        self.label_17_nmf_num.setGeometry(QtCore.QRect(25, 19, 150, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(10)
        self.label_17_nmf_num.setFont(font)
        self.label_17_nmf_num.setObjectName("label_17_nmf_num")

        self.lineEdit_8_nmf_num = QtWidgets.QLineEdit(self.frame_16)
        self.lineEdit_8_nmf_num.setGeometry(QtCore.QRect(170, 18, 30, 20))
        self.lineEdit_8_nmf_num.setPalette(palette)
        self.lineEdit_8_nmf_num.setFont(font)
        self.lineEdit_8_nmf_num.setAutoFillBackground(False)
        self.lineEdit_8_nmf_num.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.lineEdit_8_nmf_num.setText("")
        self.lineEdit_8_nmf_num.setObjectName("lineEdit_8_nmf_num")

        self.layout_14 = QtWidgets.QGridLayout(self.frame_17)

        self.label_18 = QtWidgets.QLabel(self.frame_16)
        self.label_18.setGeometry(QtCore.QRect(730, 4, 351, 41))
        sizePolicy.setHeightForWidth(self.label_18.sizePolicy().hasHeightForWidth())
        self.label_18.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.label_18.setFont(font)
        self.label_18.setObjectName("label")


        #
        #
        ###### End of NMF  #############################################

        ### Set up Labels to show the used preprocessing methods in preprocessing tab
        #
        #
        self.label_note_preprocessing_prp = QtWidgets.QLabel(self.tab1_preprocessing)
        self.label_note_preprocessing_prp.setGeometry(QtCore.QRect(5, 593, 700, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(8)
        self.label_note_preprocessing_prp.setFont(font)
        self.label_note_preprocessing_prp.setObjectName("label_note_preprocessing_prp")
        self.label_note_preprocessing_prp.setStyleSheet("background-color: rgb(255, 255, 255)")

        self.label_note_preprocessing_nmf = QtWidgets.QLabel(self.tab_6_nmf)
        self.label_note_preprocessing_nmf.setGeometry(QtCore.QRect(5, 593, 700, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(8)
        self.label_note_preprocessing_nmf.setFont(font)
        self.label_note_preprocessing_nmf.setObjectName("label_note_preprocessing_nmf")
        self.label_note_preprocessing_nmf.setStyleSheet("background-color: rgb(255, 255, 255)")


        self.label_note_preprocessing_km = QtWidgets.QLabel(self.tab_4_kmeans)
        self.label_note_preprocessing_km.setGeometry(QtCore.QRect(5, 593, 700, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(8)
        self.label_note_preprocessing_km.setFont(font)
        self.label_note_preprocessing_km.setObjectName("label_note_preprocessing_km")
        self.label_note_preprocessing_km.setStyleSheet("background-color: rgb(255, 255, 255)")
        

        self.label_note_preprocessing_pca = QtWidgets.QLabel(self.tab_2_pca)
        self.label_note_preprocessing_pca.setGeometry(QtCore.QRect(5, 593, 700, 16))
        font = QtGui.QFont()
        font.setFamily("Roboto")
        font.setPointSize(8)
        self.label_note_preprocessing_pca.setFont(font)
        self.label_note_preprocessing_pca.setObjectName("label_note_preprocessing_pca")
        self.label_note_preprocessing_pca.setStyleSheet("background-color: rgb(255, 255, 255)")
        #
        #
        #########################

        self.gridLayout_2.addWidget(self.frame_6_cls, 2, 2, 3, 1)
        self.Button_mixture = QtWidgets.QPushButton(self.frame_5_cls)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_mixture.setFont(font)
        self.Button_mixture.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_mixture.setObjectName("Button_mixture")
        self.gridLayout_2.addWidget(self.Button_mixture, 4, 0, 1, 1)
        self.Button_reference_2 = QtWidgets.QPushButton(self.frame_5_cls)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_reference_2.setFont(font)
        self.Button_reference_2.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_reference_2.setObjectName("Button_reference_2")
        self.gridLayout_2.addWidget(self.Button_reference_2, 3, 0, 1, 1)
        self.Button_run3_cls = QtWidgets.QPushButton(self.frame_5_cls)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Button_run3_cls.sizePolicy().hasHeightForWidth())
        self.Button_run3_cls.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_run3_cls.setFont(font)
        self.Button_run3_cls.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_run3_cls.setObjectName("Button_run3_cls")
        self.gridLayout_2.addWidget(self.Button_run3_cls, 2, 3, 3, 1)
        
        self.tabWidget.addTab(self.tab_6_nmf, "")
        self.tabWidget.addTab(self.tab_4_kmeans, "")
        self.tabWidget.addTab(self.tab_5_mcr, "") 
        self.tabWidget.addTab(self.tab_3_cls, "")
        
        

        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        self.tabWidget_2.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)



        ##################________________ From here the added codes for functions________________################## 


        self.Button_open.clicked.connect(self.OpenFolder)
        self.Button_save.clicked.connect(self.SaveData)
        self.Button_run1_preprocessing.clicked.connect(self.Preprocessing)
        self.Button_run2_pca.clicked.connect(self.perform_pca)
        self.checkbox_change_unit.toggled.connect(self.ChechBoxUnit)
        self.checkbox_kmean_pca.toggled.connect(self.ChechBoxUnit)
        self.checkbox_despiking.toggled.connect(self.ChechBoxUnit)

        self.Button_reference_1.clicked.connect(self.OpenFileCls)
        self.Button_reference_2.clicked.connect(self.OpenFileCls)
        self.Button_mixture.clicked.connect(self.OpenFileCls)
        self.Button_run3_cls.clicked.connect(self.perform_cls)
        self.Button_run_kmeans.clicked.connect(self.KmeansClustering)
        self.Button_kmeans_save.clicked.connect(self.SaveData)
        self.Button_run_mcr.clicked.connect(self.perform_mcr)
        self.Button_mcr_save.clicked.connect(self.SaveData)

        self.Button_run_nmf.clicked.connect(self.perform_nmf)
        self.Button_nmf_save.clicked.connect(self.SaveData)

        ###### Set up the plots of preprocessing
        self.fig1 , (self.ax1, self.ax2) = plt.subplots(2 , 1 , figsize=(10, 4), dpi=90, num = 'Raw & Preprocessed Data')
        self.fig1.subplots_adjust(bottom=0.12 ) # , hspace = 0.2
        self.ax1.set_title('Raw', y=1.0, pad=-14, fontsize=10, loc='right')
        self.ax2.set_title('Preprocessed', y=1.0, pad=-14, fontsize=10, loc='right')
        self.canvas = FigureCanvas(self.fig1)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.layout_3.addWidget(self.toolbar)
        self.layout_3.addWidget(self.canvas)

        self.Button_openseparate_preprocessing = QtWidgets.QPushButton(self.frame_3)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_openseparate_preprocessing.setFont(font)
        self.Button_openseparate_preprocessing.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_openseparate_preprocessing.setObjectName("Button_openseparate_preprocessing")
        self.Button_openseparate_preprocessing.setGeometry(QtCore.QRect(300, 14, 140, 28))
        self.Button_openseparate_preprocessing.setText("Open in Separate Windows")
        self.Button_openseparate_preprocessing.clicked.connect(self.OpenPlots_Preprocessing)

      
        ###### Set up the plots of pca
        ## Score scatter plot
        self.fig2 = plt.figure('Scatter Plot of The Score Values')
        self.ax3 = self.fig2.add_subplot(111)
        # self.ax3.set_aspect('equal', adjustable='box')
        self.ax3.set_xlabel('PC 1')
        self.ax3.set_ylabel('PC 2')
        self.canvas2 = FigureCanvas(self.fig2)
        self.toolbar2 = NavigationToolbar(self.canvas2, self)
        self.layout_4.addWidget(self.toolbar2)
        self.layout_4.addWidget(self.canvas2)

        ## Loadings
        self.fig3, (self.ax4) = plt.subplots(3, sharex=True, num = 'Mean Spectrum, 1st & 2nd Loadings' )
        self.ax4[0].set_title('Mean Spectrum', y=1.0, pad=-14, fontsize=10, loc='right')
        self.ax4[1].set_title('PC1', y=1.0, pad=-14, fontsize=10, loc='right')
        self.ax4[2].set_title('PC2', y=1.0, pad=-14, fontsize=10, loc='right')
        self.canvas3 = FigureCanvas(self.fig3)
        self.toolbar3 = NavigationToolbar(self.canvas3, self)
        self.layout_5.addWidget(self.toolbar3)
        self.layout_5.addWidget(self.canvas3)

        self.Button_openseparate_pca = QtWidgets.QPushButton(self.frame_5)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_openseparate_pca.setFont(font)
        self.Button_openseparate_pca.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_openseparate_pca.setObjectName("Button_openseparate_pca")
        self.Button_openseparate_pca.setGeometry(QtCore.QRect(300, 11, 140, 28))
        self.Button_openseparate_pca.setText("Open in Separate Windows")
        self.Button_openseparate_pca.clicked.connect(self.OpenPlots_pca)

        ## Scree plot
        self.fig4 = plt.figure('Scree Plot')
        self.ax5 = self.fig4.add_subplot(111)
        self.fig4.subplots_adjust(bottom=0.15 )
        self.ax5.set_xlabel('Principal Component')
        self.ax5.set_ylabel('Proportion of Variance Explained')
        self.canvas4 = FigureCanvas(self.fig4)
        self.toolbar4 = NavigationToolbar(self.canvas4, self)
        self.layout_6.addWidget(self.toolbar4)
        self.layout_6.addWidget(self.canvas4)

        ## Score Images
        self.fig5 , (self.ax6, self.ax7) = plt.subplots(1 , 2 , figsize=(10, 4), dpi=100, num = 'PC1 & PC2')
        self.ax6.set_title('PC 1')
        self.ax7.set_title('PC 2')
        self.canvas5 = FigureCanvas(self.fig5)
        self.toolbar5 = NavigationToolbar(self.canvas5, self)
        self.layout_7.addWidget(self.toolbar5)
        self.layout_7.addWidget(self.canvas5)

        
        ###### Set up the plots of cls
        self.fig6 = plt.figure('CLS fitting')
        self.ax8 = self.fig6.add_subplot(111)
        self.ax8.set_xlabel('Wavenumber ($cm^{-1}$)')
        self.ax8.set_ylabel('Intensity (a.u.)')        
        self.canvas6 = FigureCanvas(self.fig6)
        self.toolbar6 = NavigationToolbar(self.canvas6, self)
        self.layout_8.addWidget(self.toolbar6)
        self.layout_8.addWidget(self.canvas6)


        ###### Button for separate windows of plots of NMF
        self.Button_openseparate_nmf = QtWidgets.QPushButton(self.frame_17)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_openseparate_nmf.setFont(font)
        self.Button_openseparate_nmf.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_openseparate_nmf.setObjectName("Button_openseparate_nmf")
        self.Button_openseparate_nmf.setGeometry(QtCore.QRect(300, 11, 140, 28))
        self.Button_openseparate_nmf.setText("Open in Separate Windows")
        self.Button_openseparate_nmf.clicked.connect(self.OpenPlots_nmf)
        self.Button_openseparate_nmf.hide()

        ###### Button for separate windows of plots of K-Means
        self.Button_openseparate_kmca = QtWidgets.QPushButton(self.frame_11)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_openseparate_kmca.setFont(font)
        self.Button_openseparate_kmca.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_openseparate_kmca.setObjectName("Button_openseparate_kmca")
        self.Button_openseparate_kmca.setGeometry(QtCore.QRect(300, 11, 140, 28))
        self.Button_openseparate_kmca.setText("Open in Separate Windows")
        self.Button_openseparate_kmca.clicked.connect(self.OpenPlots_kmca)
        self.Button_openseparate_kmca.hide()

        ###### Button for separate windows of plots of mcr
        self.Button_openseparate_mcr = QtWidgets.QPushButton(self.frame_14)
        font = QtGui.QFont()
        font.setFamily("Roboto")
        self.Button_openseparate_mcr.setFont(font)
        self.Button_openseparate_mcr.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Button_openseparate_mcr.setObjectName("Button_openseparate_kmca")
        self.Button_openseparate_mcr.setGeometry(QtCore.QRect(300, 11, 140, 28))
        self.Button_openseparate_mcr.setText("Open in Separate Windows")
        self.Button_openseparate_mcr.clicked.connect(self.OpenPlots_mcr)
        self.Button_openseparate_mcr.hide()


    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "PARS"))
        self.smotthing_box.setStatusTip(_translate("MainWindow", "Smoothing via Savitzky-Golay filter"))
        self.smotthing_box.setTitle(_translate("MainWindow", "----- Smoothing ------------------------ Despiking -----"))
        self.lineEdit_4_OrderOfPoly.setText(_translate("MainWindow", "2"))
        self.lineEdit_6_Kmeans_pca.setText(_translate("MainWindow", "2"))
        self.label_9_order_of_poly.setText(_translate("MainWindow", "Order of Polynomial:"))
        self.lineEdit_3_WindowLength.setText(_translate("MainWindow", "5"))
        self.lineEdit_9_wl_despiking.setText(_translate("MainWindow", "15"))
        self.label_10_window_length.setText(_translate("MainWindow", "Window Length:"))
        self.label_20_window_length.setText(_translate("MainWindow", "Threshold:"))
        self.Button_save.setStatusTip(_translate("MainWindow", "Choose where you want to save the preprocessed data"))
        self.Button_save.setText(_translate("MainWindow", "Save"))
        self.spectral_range_box.setTitle(_translate("MainWindow", "Spectral Range"))
        self.label_5_cm.setText(_translate("MainWindow", u"cm \u207B\xb9"))
        self.label_6_cm.setText(_translate("MainWindow", u"cm \u207B\xb9"))
        self.label_7_to.setText(_translate("MainWindow", "to"))
        self.checkbox_change_unit.setText(_translate("MainWindow", u"Convert the Unit to cm \u207B\xb9"))
        self.checkbox_kmean_pca.setText(_translate("MainWindow", "K-Means Clustering of the PCA Score Values"))
        self.checkbox_despiking.setText(_translate("MainWindow", "Spike Removal"))
        self.label_8_from.setText(_translate("MainWindow", "From"))
        self.label_11_laser.setText(_translate("MainWindow", "Laser:"))
        self.label_12_nm.setText(_translate("MainWindow", "nm"))
        self.label_13_pca_clusters.setText(_translate("MainWindow", "Number of Clusters:"))
        self.Button_open.setStatusTip(_translate("MainWindow", "Choose the folder of the raw data"))
        self.Button_open.setText(_translate("MainWindow", "Open"))
        self.Button_run1_preprocessing.setStatusTip(_translate("MainWindow", "Run Preprocessing"))
        self.Button_run1_preprocessing.setText(_translate("MainWindow", "Run"))
        self.label_baseline.setText(_translate("MainWindow", "Baseline Correction Method :"))
        self.radioButton_baseline_arpls.setText(_translate("MainWindow", "arPLS"))
        self.radioButton_baseline_drpls.setText(_translate("MainWindow", "drPLS"))
        self.radioButton_baseline_none.setText(_translate("MainWindow", "None"))
        self.label_normalization.setText(_translate("MainWindow", "Normalization Method :"))
        self.radioButton_normalization_snv.setStatusTip(_translate("MainWindow", "SNV = Standard Normal Variation"))
        self.radioButton_normalization_snv.setText(_translate("MainWindow", "SNV"))
        self.radioButton_normalization_area.setText(_translate("MainWindow", "Area"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab1_preprocessing), _translate("MainWindow", "Preprocessing"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_3_scatter), _translate("MainWindow", "Scatter Plot of Score Values"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_4_loading), _translate("MainWindow", "Loadings"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_5_scree), _translate("MainWindow", "Scree Plot"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_6_pca_image), _translate("MainWindow", "PC Images"))
        self.tabWidget_kmeans.setTabText(self.tabWidget_kmeans.indexOf(self.tab_7_kmeans_image), _translate("MainWindow", "K-Means Reconstructed Image"))
        self.tabWidget_kmeans.setTabText(self.tabWidget_kmeans.indexOf(self.tab_8_kmeans_mean), _translate("MainWindow", "Mean Spectra of Clusters"))
        self.tabWidget_kmeans.setTabText(self.tabWidget_kmeans.indexOf(self.tab_9_elbow), _translate("MainWindow", "Elbow Plot"))
        self.tabWidget_mcr.setTabText(self.tabWidget_mcr.indexOf(self.tab_9_mcr_spectra), _translate("MainWindow", "MCR-ALS Retrieved Raman Spectra"))
        self.tabWidget_mcr.setTabText(self.tabWidget_mcr.indexOf(self.tab_10_mcr_image), _translate("MainWindow", "MCR-ALS Distribution Images"))
        self.label_14_kmeans_num.setText(_translate("MainWindow", "Number of Clusters:"))
        self.Button_run2_pca.setStatusTip(_translate("MainWindow", "Run Principal Component Analysis"))
        self.Button_run2_pca.setText(_translate("MainWindow", "Run"))
        self.label.setText(_translate("MainWindow", "* Before perfoming Principal Component Analysis,\n\
   preprocessing of the raw data must be done via\n   SNV normalization method."))
        self.label_15.setText(_translate("MainWindow", "* Before perfoming K-Means Cluster Analysis,\n\
   preprocessing of the raw data must be done via\n   Area normalization method."))
        self.label_16.setText(_translate("MainWindow", "* Before perfoming Multivariate Curve Resolution,\n"
                                                       "  preprocessing of the raw data via Area normalization\n"
                                                       "  method, and KMCA or SVD-NMF must be done." ))

        self.label_18.setText(_translate("MainWindow", "* Before perfoming Non-negative Matrix Factorization,\n"
                                                       "   preprocessing of the raw data must be done." ))
                                                       
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2_pca), _translate("MainWindow", "PCA"))
        self.tab_2_pca.setStatusTip(_translate("MainWindow", "PCA = Principal Component Analysis"))
        self.tab_3_cls.setStatusTip(_translate("MainWindow", "CLS = Classical Least Squares"))
        self.tab_4_kmeans.setStatusTip(_translate("MainWindow", "KMCA = K-Means Cluster Analysis"))
        self.Button_reference_1.setStatusTip(_translate("MainWindow", "Choose the first reference spectrum"))
        self.Button_reference_1.setText(_translate("MainWindow", "Reference 1"))
        self.label_2_ratio1.setText(_translate("MainWindow", "Ratio of Component 1 (%)"))
        self.label_3_ratio2.setText(_translate("MainWindow", "Ratio of Component 2 (%)"))
        self.Button_mixture.setStatusTip(_translate("MainWindow", "Choose the spectrum of the mixture sample"))
        self.Button_mixture.setText(_translate("MainWindow", "Mixture"))
        self.Button_reference_2.setStatusTip(_translate("MainWindow", "Choose the second reference spectrum"))
        self.Button_reference_2.setText(_translate("MainWindow", "Reference 2"))
        self.Button_run3_cls.setStatusTip(_translate("MainWindow", "Run Classical Least Squares Method"))
        self.Button_run3_cls.setText(_translate("MainWindow", "Run"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3_cls), _translate("MainWindow", "CLS"))
        self.Button_run_kmeans.setText(_translate("MainWindow", "Run"))
        self.Button_run_kmeans.setStatusTip(_translate("MainWindow", "Run K-Means Cluster Analysis"))
        self.Button_kmeans_save.setText(_translate("MainWindow", "Export Mean\n Spectra"))
        self.Button_kmeans_save.setStatusTip(_translate("MainWindow", "Save the mean spectra of the clusters in separate text files"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_5_mcr), _translate("MainWindow", "MCR-ALS"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_6_nmf), _translate("MainWindow", "SVD-NMF"))
        self.Button_run_mcr.setText(_translate("MainWindow", "Run"))
        self.Button_run_mcr.setStatusTip(_translate("MainWindow", "Run Multivariate Curve Resolution - Alternating Regression"))
        self.Button_mcr_save.setText(_translate("MainWindow", "Export \nSpectra"))
        self.Button_mcr_save.setStatusTip(_translate("MainWindow", "Save the spectra of pure components in separate text files"))
        self.Button_nmf_save.setText(_translate("MainWindow", "Export NMF\nComponents"))
        self.Button_nmf_save.setStatusTip(_translate("MainWindow", "Save the decomposed spectra in separate text files"))
        self.Button_run_nmf.setText(_translate("MainWindow", "Run"))
        self.Button_run_nmf.setStatusTip(_translate("MainWindow", "Run Non-negative Matrix Factorization"))
        self.label_17_nmf_num.setText(_translate("MainWindow", "Number of Components:"))
        self.lineEdit_8_nmf_num.setText(_translate("MainWindow", ""))
        self.label_19.setText(_translate("MainWindow", "Initial Spectra:"))
        self.radioButton_initial_nmf.setText(_translate("MainWindow", "SVD-NMF"))
        self.radioButton_initial_kmeans.setText(_translate("MainWindow", "KMCA"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_4_kmeans), _translate("MainWindow", "KMCA"))

        ##### Added status tips #####
        self.lineEdit_1_range.setStatusTip(_translate("MainWindow", "If you leave it empty, it is considered 0"))
        self.lineEdit_5_laser.setStatusTip(_translate("MainWindow","Excitation Wavelength"))
        self.lineEdit_2_range.setStatusTip(_translate("MainWindow", "If you leave it empty, it is considered until the end"))
        self.lineEdit_3_WindowLength.setStatusTip(_translate("MainWindow", "Window Length = Number of Points,\
 The length of the window for smoothing must be an odd number"))
        self.label_10_window_length.setStatusTip(_translate("MainWindow", "Window Length = Number of Points,\
 The length of the window for smoothing must be an odd number"))
        self.label_9_order_of_poly.setStatusTip(_translate("MainWindow", "Smoothing via Savitzky-Golay filter"))
        self.label.setStatusTip(_translate("MainWindow", "PCA does not work if preprocessing is not performed yet"))
        self.radioButton_baseline_drpls.setStatusTip(_translate("MainWindow", "drPLS  = Doubly Reweighted Penalized Least Squares"))
        self.radioButton_baseline_arpls.setStatusTip(_translate("MainWindow", "arPLS = Asymmetrically Reweighted Penalized Least Squares Smoothing"))
        self.radioButton_baseline_none.setStatusTip(_translate("MainWindow", "None = No Baseline Correction"))
        self.frame_6_cls.setStatusTip(_translate("MainWindow", "The automatic preprocessing in this part includes zero correction,\
 area normalization and baseline correction via arPLS method"))
        self.radioButton_pcacolormap.setStatusTip(_translate("MainWindow", "Unchecked = jet , Checked = Greys"))       
        self.checkbox_despiking.setStatusTip(_translate("MainWindow", "Removes spikes in spectra that made by Cosmic rays (Muons)"))
        

    def OpenFolder(self): #imports the raw data and puts them in a hypercube

        msg = QMessageBox()
        items = ["Txt (Single File)","Txt (Image)", "HDF5 (Image)" ]
        self.data_format, okpressed = QtWidgets.QInputDialog.getItem(self, 'Format', 'Format of the Raw Data:' , items\
            , 0 , False)
        QApplication.quitOnLastWindowClosed() == False

        if self.data_format == 'Txt (Single File)':
            txt_file_path= QFileDialog.getOpenFileName(self, 'Select the Text File', filter = "Text files (*.txt *.asc)")
            if txt_file_path[0] == '':
                pass
            
            else:
                df= pd.read_table(txt_file_path[0], header=0 , decimal="," , index_col=False , skiprows=34)
                self.x_ref1 = df.iloc[:,0]
                self.x_ref1 = self.x_ref1.values
                self.x_ref1 = self.x_ref1.astype(np.float)
                y_ref1 = df.iloc[:,1]
                y_ref1 = y_ref1.values
                y_ref1 = y_ref1.astype(np.float)
                self.size_x = 1
                self.size_y = 1
                self.hypercube = np.zeros((len(self.x_ref1),self.size_x,self.size_y))
                self.hypercube[: , 0 , 0] = y_ref1
                self.check_open = 1 

                # with open('txt_file_path') as f:
                #     lines=f.readlines()
                #     result=[]
                #     for x in lines:
                #         result.append(x.split(' ')[0])


        elif self.data_format == 'Txt (Image)':
            directory_path_open = str(QFileDialog.getExistingDirectory(self, "Select the Directory of the Txt Files"))
            if directory_path_open == '':
                pass
            else:
                list_of_files = os.listdir(directory_path_open) 
                
                ##### Hypercube #####

                ## 1) Find the size of the cube via checking the names of the data

                x_00_indices_list=[]
                y_00_indices_list=[]


                x00_axis_list = [file for file in list_of_files if "[00 x" in file]
                y00_axis_list = [file for file in list_of_files if "x 00]" in file]

                for name in x00_axis_list:
                    for idx, j in enumerate(name):
                        if name[idx:idx+5] == '[00 x':
                            y_00_indices_list.append(int(name[idx+6:idx+8])+1)

                for name in y00_axis_list:
                    for idx, j in enumerate(name):
                        if name[idx:idx+5] == 'x 00]':
                            x_00_indices_list.append(int(name[idx-3:idx-1])+1)


                ### If the name of the single file in the folder does not include [aa x bb] , For opening one spectrum
                if len(x_00_indices_list) == 0 or len(y_00_indices_list) == 0:

                    for file in list_of_files:
                        if fnmatch.fnmatch(file, '*.asc') or fnmatch.fnmatch(file, '*.txt'):
                            df2= pd.read_table(directory_path_open + '/' + file, header=0 , decimal="," , index_col=False , skiprows=34)
                            self.x_ref1 = df2.iloc[:,0]
                            self.x_ref1 = self.x_ref1.values
                            y_ref1 = df2.iloc[:,1]
                            y_ref1 = y_ref1.values
                            y_ref1 = y_ref1.astype(np.float)
                            self.hypercube = np.zeros((len(self.x_ref1), 1 ,1))
                            self.hypercube[:,0,0] = y_ref1
                            self.size_x = 1
                            self.size_y = 1
                            pass

                else:
                    if len(x_00_indices_list) == 1: # For the case there is one mosaic [00 x 00]
                        self.size_x = 1
                    elif len(x_00_indices_list) >= 1:
                        self.size_x = max(x_00_indices_list) # Number of the rows

                    if len(y_00_indices_list) == 0: # For the case there is one mosaic [00 x 00]
                        self.size_y = 1
                    elif len(y_00_indices_list) >= 1:
                        self.size_y = max(y_00_indices_list) # Number of the columns


                    for file in list_of_files:
                        if fnmatch.fnmatch(file, '*.asc') or fnmatch.fnmatch(file, '*.txt'):
                            df3= pd.read_table(directory_path_open + '/' + file, header=0 , decimal="," , index_col=False , skiprows=34)
                            self.x_ref1 = df3.iloc[:,0]
                            self.x_ref1 = self.x_ref1.values
                            self.x_ref1 = self.x_ref1.astype(np.float)
                            y_ref1 = df3.iloc[:,1]
                            y_ref1 = y_ref1.values
                            y_ref1 = y_ref1.astype(np.float)
                            break


                    ## 2) Make a Zeros 3d matrix (Cube) equal to the size of the data (Number of mosaics)
                    self.hypercube = np.zeros((len(self.x_ref1),self.size_x,self.size_y)) 

                    ## 3) Put the data in the right position in the cube according to the names of the data

                    for file in list_of_files:
                        if fnmatch.fnmatch(file, '*.asc') or fnmatch.fnmatch(file, '*.txt'):
                            df4= pd.read_table(directory_path_open + '/' + file, header=0 , decimal="," , index_col=False , skiprows=34)
                            self.x_ref1 = df4.iloc[:,0]
                            self.x_ref1 = self.x_ref1.values
                            self.x_ref1 = self.x_ref1.astype(np.float)
                            y_ref1 = df4.iloc[:,1]
                            y_ref1 = y_ref1.values
                            y_ref1 = y_ref1.astype(np.float)
                            for idx, letter in enumerate(file):
                                for a in range(0,10):
                                    for b in range(0,10):
                                        if file[idx:idx+5] == f'[{a}{b} x':
                                            index_x = int(f'{a}{b}')
                                            break 
                                    #break

                                for c in range(0,10):
                                    for d in range(0,10):
                                        if file[idx:idx+5] == f'x {c}{d}]':
                                            index_y = int(f'{c}{d}')
                                            print(index_y)
                                            break
                                    #break   

                            # print(index_x , index_y)
                            self.hypercube[:,index_x, index_y] = y_ref1
                            # print(hypercube[:,index_x, index_y])

                    self.check_open = 1                 
                    print(list_of_files)

        elif self.data_format == 'HDF5 (Image)':

            hdf5_file_path= QFileDialog.getOpenFileName(self, 'Select the HDF5 File')
            if hdf5_file_path[0] == '':
                pass
            
            else:
                hdf = h5py.File(hdf5_file_path[0] , 'r')
                data1 = hdf.get('hypercube1')
                self.hypercube = np.array(data1)

                data2= hdf.get('wavenumber')
                self.x_ref1 = np.array(data2)
                hdf.close()

                self.size_x = len(self.hypercube[0])
                self.size_y = len(self.hypercube[0][0])
                self.check_open = 1    

    def ChechBoxUnit (self): # To activate the laser wavelength input for unit change, or kmeans clusterig of PCA results, or despiking
        if self.checkbox_change_unit.isChecked() == True:
            self.label_12_nm.setEnabled(True)
            self.lineEdit_5_laser.setEnabled(True)
            self.label_11_laser.setEnabled(True)
        else:
            self.label_12_nm.setEnabled(False)
            self.lineEdit_5_laser.setEnabled(False)
            self.label_11_laser.setEnabled(False)

        if self.checkbox_kmean_pca.isChecked() == True:
            self.lineEdit_6_Kmeans_pca.setEnabled(True)
            self.label_13_pca_clusters.setEnabled(True)
        else:
            self.lineEdit_6_Kmeans_pca.setEnabled(False)
            self.label_13_pca_clusters.setEnabled(False)

        if self.checkbox_despiking.isChecked() == True:
            self.lineEdit_9_wl_despiking.setEnabled(True)
            self.label_20_window_length.setEnabled(True)
        else:
            self.lineEdit_9_wl_despiking.setEnabled(False)
            self.label_20_window_length.setEnabled(False)    


    def UnitChange(self):      
                                       
        if self.checkbox_change_unit.isChecked() == True:
            if self.laser_wavelength == '':
                self.laser_wavelength, okpressed = QtWidgets.QInputDialog.getText(self, 'Error', 'Please enter the wavelength of the laser (nm):')
                self.lineEdit_5_laser.setText(self.laser_wavelength)
                if self.laser_wavelength == '': # For the case either cancel is clicked or nothing is entered
                    msg = QMessageBox()
                    msg.setWindowTitle("FYI")
                    msg.setText(" -> Wavelength = 777 nm :) ")
                    msg.setIcon(QMessageBox.Information) 
                    msg.setStandardButtons(QMessageBox.Ok)
                    x = msg.exec_()
                    if x == QMessageBox.Ok:
                        self.laser_wavelength = '777'                       
            self.x_ref1_wn= 10**7 * ((1/int(self.laser_wavelength))-(1/self.x_ref1))
            self.note_unit_conversion = 'Unit Conversion: Yes'
        else:
            self.note_unit_conversion = 'Unit Conversion: No'
            self.x_ref1_wn= self.x_ref1
            pass
            

        self.check_open == 0
        return(self.x_ref1_wn)
         

    def spectral_range (self):

        self.roi_1 = self.lineEdit_1_range.text()
        self.roi_2 = self.lineEdit_2_range.text()

        if self.roi_1 == '':
            self.roi_1 = 0
        if self.roi_2 == '':
            self.roi_2 = self.x_ref1_wn[-1]

        if int(self.roi_2) - int(self.roi_1) < 14:
            msg = QMessageBox()
            msg.setWindowTitle("Error")
            msg.setText("The spectral range must be larger than 14 cm \u207B\xb9 !")
            msg.setInformativeText("The end point is added with 14.")
            msg.setIcon(QMessageBox.Critical) 
            msg.setStandardButtons(QMessageBox.Ok)
            x = msg.exec_()
            if x == QMessageBox.Ok:
                self.roi_2 = int(self.roi_2) + 14
                self.lineEdit_2_range.setText(f'{int(self.roi_2)}')

        y_index_list=[]
        self.y_ref1_new= np.array([])
        self.x_ref1_new = self.x_ref1_wn[ (int(self.roi_1) <= self.x_ref1_wn) & (self.x_ref1_wn <= int(self.roi_2))]

        for i in self.x_ref1_new:
            y_index = np.where(self.x_ref1_wn == i)
            y_index_list.append(int(y_index[0]))

        for j in y_index_list:
            self.y_ref1_new = np.append(self.y_ref1_new, self.spectrum1[j])

        return (self.x_ref1_new,self.y_ref1_new)


    def correct_zero (self, y_data):

        y_ref1_min = min(y_data) 
        if y_ref1_min < 0:
            for i in range (0,len(y_data)):
                y_data[i] = y_data[i] + abs(y_ref1_min)
        else:
            for i in range (0,len(y_data)):
                y_data[i] = y_data[i] - abs(y_ref1_min)

        return (y_data)


    def despiking (self , y_data1):


        #____Whitaker and Hayes’:

        threshold_ = int(self.lineEdit_9_wl_despiking.text())

        def modified_z_score(intensity):
            median_int = np.median(intensity)
            mad_int = np.median([np.abs(intensity - median_int)])
            modified_z_scores = 0.6745 * (intensity - median_int) / mad_int
            return modified_z_scores    

        dist = 0
        delta_intensity = [] 
        delta_int = np.diff(y_data1)
        # Alternatively to the for loop one can use: 
        # delta_int = np.diff(intensity)

        threshold = threshold_


        intensity_modified_z_score = np.array(np.abs(modified_z_score(delta_int)))

        spikes = abs(intensity_modified_z_score) > threshold

        # for i in intensity_modified_z_score:
        #         if abs(intensity_modified_z_score) > threshold:
        #                 spikes = intensity_modified_z_score.index(i)
        y_out = y_data1.copy()

        for i in np.arange(len(spikes)):
                if spikes[i] != 0: # If we have an spike in position i
                        w = np.arange(i-10,i+1+10) # we select 2 m + 1 points around our spike
                        w2 = w[spikes[w] == 0] # From such interval, we choose the ones which are not spikes
                        y_out[i] = np.mean(y_data1[w2]) # and we average their values


        ###### Second Despiker (If necessary it could be added to the first one. Its speed is lower)

        # size = int(self.lineEdit_9_wl_despiking.text())
        # despiked_y_data = np.array([])


        # def box(A, j, size=1, nan=False):

        #     nx = A.shape[0]
        #     t = j - size if j - size > 0 else 0
        #     b = j + size if j + size < nx else nx

        #     return A[t:b+1]


        # def mean(A):
        #     return np.nanmean(A), np.nanstd(A)


        # def median(A):
        #     return np.nanmedian(A), np.nanstd(A)


        # def spikes(y_out, method='mean', size=6, n=3):
        #     '''
        #     Search spikes in the image using a moving box 
        #     of size `size` using `method` method with ± n × std
        #     '''

        #     _spikes = np.zeros((len(y_out)))
        #     for j in range(len(y_out)):
        #         m, s = mean(box(y_out, j, size, nan=True))
        #         if y_out[j] < m - n * s or y_out[j] > m + n * s:
        #             _spikes[j] = 1
        #             # print(j)
                        
        #     return _spikes == 1


        # def fill(y_out, mask, method='mean', size = 6):
        #     '''
        #     Fill the masked pixels with the mean/median value of
        #     surrounding neighboors inside a box of size `size`
        #     '''
        #     A = np.copy(y_out)
        #     nx = A.shape[0]
        #     x = np.array(range(nx))

        #     for x in x[mask]:
        #         A[x] = np.nanmedian(box(A, x , size=size, nan=True))

        #     return A


        # def clean(y_out, mask='mean', size=2, n=3, fill_method='median', fill_size=10):
        #     '''
        #     Clean image from spikes with the `fill_method` method with surrounding
        #     neighboors inside a box of size `fill_size`.
        #     '''
        #     return fill(y_out, spikes(y_out, mask, size, n), fill_method, fill_size)
        
        # spikes(y_out , size = size ) # Search the location of spikes in the image
        # despiked_y_data = clean(y_out , size = size ) # Clean the image from spikes

        
                        
        return(y_out)


    def snv_normalization (self):

        y_bc_mean= st.mean(self.y_baseline_corrected1.flatten()) #bc: Baseline Corrected
        y_bc_stdev= st.stdev(self.y_baseline_corrected1.flatten())
        self.y = (self.y_baseline_corrected1 - y_bc_mean)/y_bc_stdev

        return (self.y)


    def color_maker (self): # Makes a color for e.g. plots
        
        r= random.randrange(0, 100)
        r = r/100

        g= random.randrange(0, 100)
        g = g/100

        b= random.randrange(0, 100)
        b = b/100

        self.color = (r,g,b)

        return (self.color)

    def SaveData (self):

        sender = self.sender()
        if sender.text() == 'Save':
            self.directory_path_save = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
            if self.directory_path_save == '':
                self.check_save_preprocessing = 0
                pass
            else:
                self.check_save_preprocessing = 1
                return(self.directory_path_save)

        elif sender.text() == 'Export Mean\n Spectra':
            self.directory_path_save_kmean = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
            if self.directory_path_save_kmean == '':
                pass
            else:
                self.check_save_kmeans = 1
                return(self.directory_path_save_kmean)
            pass

        elif sender.text() == 'Export \nSpectra':
            self.directory_path_save_mcr = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
            if self.directory_path_save_mcr == '':
                pass
            else:
                self.check_save_mcr = 1
                return(self.directory_path_save_mcr)
            pass

        elif sender.text() == 'Export NMF\nComponents':
            self.directory_path_save_nmf = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
            if self.directory_path_save_nmf == '':
                pass
            else:
                self.check_save_nmf = 1
                return(self.directory_path_save_nmf)
            pass
    

    ######## Functions for openning separate windows for plots:
    #
    #
    #

    def OpenPlots_Preprocessing(self):

        if self.check_run_preprocessing == 0:
            pass
        else:

            fig_raw = plt.figure('Raw' , figsize=(12, 5))
            for counter1 in range(0,self.size_x):
                    for counter2 in range(0,self.size_y):           
                        plt.plot(self.x_ref1 , self.hypercube[:,counter1,counter2] )
            
            fig_pre = plt.figure('Preprocessed' , figsize=(12, 5))
            for i in range(0,len(self.pdm)):           
                plt.plot(self.x_ref1_new , self.pdm[i,:])
                plt.xlabel('$cm^{-1}$', fontsize=10 )

            fig_pre.show()
            fig_raw.show()

    def OpenPlots_pca(self):

        if self.check_pca == 0:
            pass
        else:         

            fig_mean = plt.figure('Mean Spectrum' , figsize=(12, 5))
            plt.plot(self.x_ref1_new, self.pdm_means, 'tab:blue' )
            
            fig_pca1 = plt.figure('PC1' , figsize=(12, 5))
            plt.plot(self.x_ref1_new, self.loadings[:,0], 'tab:red')

            fig_pca2 = plt.figure('PC2' , figsize=(12, 5))
            plt.plot(self.x_ref1_new, self.loadings[:,1], 'tab:green')

            fig_mean.show()
            fig_pca1.show()
            fig_pca2.show()


    def OpenPlots_nmf(self):

        if self.check_nmf == 0:
            pass
        else:         

            for j in range(0, int(self.num_nmf_comp)):

                fig_nmf = plt.figure(f'NMF {j+1}' , figsize=(12, 5))
                plt.plot(self.x_ref1_new, self.nmf_comp[j,:] )
                fig_nmf.show()

    def OpenPlots_kmca(self):

        if self.check_kmca == 0:
            pass
        else:   
            
            for k in range(0,np.max(self.kmeans_labels)+1):
                print(k)
                for j in range(0, len(self.kmeans_labels)):
                    if self.kmeans_labels[j] == k :
                        print(self.kmeans_labels[j])
                        fig_kmca = plt.figure(f'KMCA {self.kmeans_labels[j]+1}' , figsize=(12, 5))   
                        plt.plot(self.x_ref1_new, self.kmeans_cluster_centers[self.kmeans_labels[j],:])
                        fig_kmca.show()
                        break                   
                    else:
                        continue


    def OpenPlots_mcr(self):

        if self.check_mcr_clear == 0 :
            pass
        else:         
            
            #### Plot Spectra , Save
            for j in range(self.mcrar_shape[1]):
                fig_mcr = plt.figure(f'MCR {j+1}' , figsize=(12, 5))
                plt.plot(self.x_mcr, self.mcrar_ST_opt[j,:])
                fig_mcr.show()


    #
    #
    ##################################


    def Preprocessing(self): #connected to run
        self.laser_wavelength = self.lineEdit_5_laser.text()
        self.radioButton_normalization_snv.setStyleSheet("")
        self.radioButton_normalization_area.setStyleSheet("")

        # Check normalization for pca
        if self.radioButton_normalization_area.isChecked() == True:
            self.norm_area_check = 1
        elif self.radioButton_normalization_area.isChecked() == False:
            self.norm_area_check = 0

        if self.check_open == 1:

            self.final_data_list = []
            self.ax1.clear()
            self.ax2.clear()
            

            wl = self.lineEdit_3_WindowLength.text()
            oop = self.lineEdit_4_OrderOfPoly.text()

            if wl == '': # If the line edits are left empty
                wl = 3
            if oop == '':
                oop = 2

            if int(self.lineEdit_9_wl_despiking.text()) <= 0:
                self.lineEdit_9_wl_despiking.setText("1")


            if (int(wl) % 2) == 0 or int(oop) >= int(wl): # Check if window length is odd
                msg1 = QMessageBox()
                msg1.setWindowTitle("Error")
                msg1.setText("The length of the window for smoothing must be an odd number.\nAnd the window length must be larger"\
                    " than the order of the polynomial.")
                msg1.setIcon(QMessageBox.Critical) 
                msg1_b = msg1.setStandardButtons(QMessageBox.Ok)
                return_value = msg1.exec_()
                if return_value == QMessageBox.Ok:
                    pass
            else:
                ##### Preprocessing #####
                for counter1 in range(0,self.size_x):
                    for counter2 in range(0,self.size_y):

                        self.spectrum = np.array(self.hypercube[:,counter1,counter2])

                        #### 1) Change the unit to wavenumber

                        self.x_ref1_wn = self.UnitChange()


                        #### 2) Cosmic Ray Spike Removal

                        self.spectrum1 = self.spectrum.copy()

                        if self.checkbox_despiking.isChecked() == True:
                            self.spectrum1[:] = self.despiking (self.spectrum[:])
                            self.note_spike_removal = 'Spike Removal: yes'
                        elif self.checkbox_despiking.isChecked() == False:
                            self.note_spike_removal = 'Spike Removal: No'

                        #### 3) Truncate the Spectral Range 

                        self.x_ref1_new,self.y_ref1_new0 = self.spectral_range ()   


                        #### 4) Correct Zero of the Spectra 

                        # self.y_ref1_new = self.correct_zero (self.y_ref1_new0)
                        self.y_ref1_new = self.y_ref1_new0


                        #### 5) Smoothing 

                        y_sg_filtered1 = ss.savgol_filter(self.y_ref1_new,int(wl),int(oop)) ##Suggested in a paper: wl = 7 , oop = 5
                        self.note_smoothing = f'Smoothing: WL:{int(wl)}  OOP:{int(oop)}'

                        # y_sg_filtered1= signal.wiener(self.y_ref1_new,int(wl))

                        
                        #### 6) Baseline Correction                   

                        baseline_roi1 = []
                        for i in range (0,4000,500):
                            a = [i,i+500]
                            baseline_roi1.append(a)

                        baseline_roi = np.array(baseline_roi1)
                        #y_baseline_corrected, background = rp.baseline(x_ref1_new , y_sg_filtered1, baseline_roi, 'poly',polynomial_order=3 )

                        self.y_baseline_corrected1 = np.array([])

                        if self.radioButton_baseline_arpls.isChecked() == True:
                            y_baseline_corrected, base_arpsl = rp.baseline(self.x_ref1_new , y_sg_filtered1, \
                                baseline_roi,'arPLS',lam=10**8,ratio=0.01)
                            for i in range (len(y_baseline_corrected)):
                                self.y_baseline_corrected1 = np.append(self.y_baseline_corrected1, y_baseline_corrected[i][0])
                            self.note_baseline= 'Baseline Correction: arPLS' #For writing preprocesing specifications.


                        elif self.radioButton_baseline_drpls.isChecked() == True:
                            y_baseline_corrected, base_drpsl = rp.baseline(self.x_ref1_new , y_sg_filtered1, baseline_roi,'drPLS' ,lam=10**4)
                            for i in range (len(y_baseline_corrected)):
                                self.y_baseline_corrected1 = np.append(self.y_baseline_corrected1, y_baseline_corrected[i][0])
                            self.note_baseline= 'Baseline Correction: drPLS'
                            

                        elif self.radioButton_baseline_none.isChecked() == True:
                            self.y_baseline_corrected1 = np.array(y_sg_filtered1)
                            self.note_baseline= 'Baseline Correction: None'

                                               

                        #### 7) Normalization
                        # if self.check_pca == 0:
                        if self.radioButton_normalization_snv.isChecked() == True:
                            # ## SNV Normalization
                            # self.y_normed = self.snv_normalization()
                            # self.note_normalization= 'Normalization: SNV'
                            
                            self.y_normed = self.y_baseline_corrected1

                        elif self.radioButton_normalization_area.isChecked() == True:
                            # Area Normalization
                            self.y_normed = rp.normalise(self.y_baseline_corrected1,self.x_ref1_new,method="area")
                            self.y_normed = self.correct_zero (self.y_normed) #Correct Zero
                            self.note_normalization= 'Normalization: Area'

                            # self.y_normed = self.y_baseline_corrected1

                        # elif self.check_pca ==1:
                        #     self.y_normed = self.snv_normalization()



                        #### 8) Plot The Spectra 

                        # if self.check_pca == 1:
                        #     pass
                        # else:
                        # color = self.color_maker()                        
                        cmap = plt.cm.get_cmap('plasma')
                        b= random.randrange(0, 100)
                        b = b/100


                        # plt.figure('Raw')
                        self.ax1.plot(self.x_ref1 , self.spectrum, c=cmap(b)) #label=file #c=color
                        
                        if self.checkbox_change_unit.isChecked():
                            self.ax1.set_xlabel('nm', horizontalalignment='right', x=1.0 , y=1.0)
                            self.ax1.xaxis.set_label_coords(1.03, -0.015)
                        else:
                            self.ax1.set_xlabel('$cm^{-1}$', fontsize=10)
                            self.ax1.xaxis.set_label_coords(1.03, -0.015)


                        # plt.figure('Preprocessed')
                        self.ax2.plot(self.x_ref1_new , self.y_normed, c=cmap(b)) #label=file
                        self.ax2.set_xlabel('$cm^{-1}$' , fontsize=10)
                        self.ax2.xaxis.set_label_coords(1.03, -0.015)

                        self.ax1.set_title('Raw',y=1.0, pad=-14, fontsize=10, loc='right')
                        self.ax2.set_title('Preprocessed',y=1.0, pad=-14, fontsize=10, loc='right')

                        # plt.legend(bbox_to_anchor=(1.2, 1), loc="upper right", fontsize='x-small')
                        #plt.xlim([int(x_min), int(x_max)]) #Change the range on the plot
                        #plt.xlim([int(x_min), int(x_max)]) #Change the range on the plot
                        
                        

                        #### 9) Save The Processed Data
                        # if self.check_pca == 1:
                        #     pass
                        # else:                            
                        if self.check_save_preprocessing == 1:
                            if self.data_format == 'Txt (Single File)':
                                df = pd.DataFrame(np.column_stack((self.x_ref1_new, self.y_normed)))
                                df.to_csv(self.directory_path_save + '/' + 'New_' + '.txt',\
                                 sep='\t', header=False, index=False)

                            if self.data_format == 'Txt (Image)':
                                if counter1 <= 9 and counter2 <= 9:
                                    df = pd.DataFrame(np.column_stack((self.x_ref1_new, self.y_normed)))
                                    df.to_csv(self.directory_path_save + '/' + 'New_' + f'[0{counter1} x 0{counter2}]' + '.txt',\
                                     sep='\t', header=False, index=False)
                                elif counter1 <= 9 and counter2 > 9:
                                    df = pd.DataFrame(np.column_stack((self.x_ref1_new, self.y_normed)))
                                    df.to_csv(self.directory_path_save + '/' + 'New_' + f'[{counter1} x 0{counter2}]' + '.txt',\
                                     sep='\t', header=False, index=False)
                                elif counter1 > 9 and counter2 <= 9:
                                    df = pd.DataFrame(np.column_stack((self.x_ref1_new, self.y_normed)))
                                    df.to_csv(self.directory_path_save + '/' + 'New_' + f'[0{counter1} x {counter2}]' + '.txt',\
                                     sep='\t', header=False, index=False)
                                elif counter1 > 9 and counter2 > 9:
                                    df = pd.DataFrame(np.column_stack((self.x_ref1_new, self.y_normed)))
                                    df.to_csv(self.directory_path_save + '/' + 'New_' + f'[{counter1} x {counter2}]' + '.txt',\
                                     sep='\t', header=False, index=False)

                            elif self.data_format == 'HDF5 (Image)' :
                                pass

                        else:
                            pass


                        self.final_data_list.append(self.y_normed)
                        print('Done!')

                #### 10) Make a new hypercube, store the new data in it, and save it as HDF5 file

                self.pdm = np.array(self.final_data_list) # pdm = preprocessed data matrix: Prepares the data matrix for the further data analysis
                hypercube_new = np.zeros((len(self.y_normed), self.size_x , self.size_y)) # Hypercube Preprocessed
                
                ### Convert pdm to hypercube_new to save as HDF5
                i = 0
                for counter3 in range(0,self.size_x):
                    for counter4 in range(0,self.size_y):
                        hypercube_new[:, counter3 , counter4] = self.pdm[i, :]
                        i += 1

                self.pdm = np.nan_to_num(self.pdm)

                ####Save for Matlab
                # scipy.io.savemat(self.directory_path_save + '/' + 'out.mat', {'mydata': self.pdm})

                ### Save the HDF5 File
                if self.check_save_preprocessing == 1:
                    if self.data_format == 'HDF5 (Image)':
                        hdf = h5py.File(self.directory_path_save + '/' + 'Preprocessed_Data.h5' , 'w')
                        hdf.create_dataset('wavenumber' , data = self.x_ref1_new )
                        hdf.create_dataset('hypercube1' ,  data = hypercube_new)


                #### 11) show the plots on the canvas + preprocessing info

                # if self.check_pca == 0:
                self.canvas.draw()
                # plt.show()
                self.check_run_preprocessing = 1
                self.label_note_preprocessing_prp.setText(
                 f"Preprocessing: {self.note_unit_conversion} | {self.note_spike_removal} | {self.note_smoothing} | "
                 f"{self.note_baseline} | {self.note_normalization}")

                return(self.final_data_list)

                
        else:
            pass
        

    def perform_pca(self):

        if self.check_run_preprocessing == 0:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information) 
            msg.setWindowTitle("Note")
            msg.setText("Please do preprocessing (with SNV normalization), before Performing Principle Component Analysis (PCA).")
            msg.setStandardButtons(QMessageBox.Ok)
            x = msg.exec_()
            if x == QMessageBox.Ok:
                pass
        else:

            if self.size_x == 1 and self.size_y == 1:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Information) 
                msg.setWindowTitle("Note")
                msg.setText("For performing PCA more than one Raman data is needed.")
                msg.setStandardButtons(QMessageBox.Ok)
                x = msg.exec_()
                if x == QMessageBox.Ok:
                    pass

            else:
                # self.check_pca = 1
                self.ax3.clear()
                for i in range(0,3):
                    self.ax4[i].clear() 
                self.ax5.clear()
                self.ax6.clear()
                self.ax7.clear()
                if self.check_colorbar == 1:
                    self.cb1.remove()
                    self.cb2.remove()
                    self.check_colorbar == 0
                

                # self.final_data_list = self.Preprocessing() 
                # To be decided -> without this line the speed of execution
                # of PCA gets much more but with this line there is no need for running preprocessing always before PCA

                if self.norm_area_check == 1: # if preprocessing is done with area normalization
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Information) 
                    msg.setWindowTitle("Important")
                    msg.setText("Please note that for PCA, the normalization of raw data must be done via SNV method!")
                    msg.setInformativeText("If you performed preprocessing with area normalization, the PCA results will not be valid!"\
                                            " Please repeat the preprocessing with SNV normalization and then run PCA again.")
                    self.radioButton_normalization_area.setChecked(False)
                    self.radioButton_normalization_snv.setChecked(True)
                    self.radioButton_normalization_snv.setStyleSheet("QRadioButton::indicator:checked"
                                        "{background-color:blue ; border:                 2px solid white;}"
                                        "QRadioButton" 
                                        "{background-color:240,240,240 ; color:                  black;}"
                                        "QRadioButton::indicator" 
                                        "{width:                  10px; height:                 10px; border-radius:          7px;}"
                                        "QRadioButton::indicator:unchecked {background-color:       white; border: 2px solid white;}"
                                        ) 
                    msg.setStandardButtons(QMessageBox.Ok)
                    x = msg.exec_()
                    if x == QMessageBox.Ok:
                        pass

                # pdm = np.array(self.final_data_list) # pdm = preprocessed data matrix: Prepares the data matrix for the further data analysis

                pca = decomposition.PCA() # n_components = 2
                pca.fit(self.pdm)
                scores = pca.transform(self.pdm)
                self.loadings = pca.components_.T * np.sqrt(pca.explained_variance_)

                ##### Score value scatter plot
                pca_cluster_num = self.lineEdit_6_Kmeans_pca.text()

                # Color coding the scatter plot via k-means
                if self.checkbox_kmean_pca.isChecked() == True and pca_cluster_num != '':
                    kmeans_pca = KMeans(n_clusters = int(pca_cluster_num))
                    kmeans_pca.fit(scores)
                    self.ax3.scatter(scores[:,0],scores[:,1], c=kmeans_pca.labels_, cmap='rainbow')
                    self.ax3.set_xlabel('PC 1')
                    self.ax3.set_ylabel('PC 2')
                    self.canvas2.draw()

                elif self.checkbox_kmean_pca.isChecked() == False or pca_cluster_num == '': 
                    self.ax3.scatter(scores[:,0],scores[:,1])
                    self.ax3.set_xlabel('PC 1')
                    self.ax3.set_ylabel('PC 2')
                    self.canvas2.draw()

                
                # for i, txt in enumerate(list_of_files):
                #   plt.annotate(f'{txt}', (scores[i,0] ,scores[i,1] ))

                ###### Loadings plots
                self.pdm_means = self.pdm.mean(axis=0)
                
                self.ax4[0].plot(self.x_ref1_new, self.pdm_means, 'tab:blue')
                self.ax4[0].set_title('Mean Spectrum', y=1.0, pad=-14, fontsize=10, loc='right')

                self.ax4[1].plot(self.x_ref1_new, self.loadings[:,0], 'tab:red')
                self.ax4[1].set_title('PC1', y=1.0, pad=-14, fontsize=10, loc='right')

                self.ax4[2].plot(self.x_ref1_new, self.loadings[:,1], 'tab:green')
                self.ax4[2].set_title('PC2', y=1.0, pad=-14, fontsize=10, loc='right')
                self.ax4[2].set(xlabel='Raman Shift ($cm^{-1}$)')
                self.canvas3.draw()
                
                ##### Scree plot
                PC_values = np.arange(pca.n_components_) + 1
                self.ax5.plot(PC_values, pca.explained_variance_ratio_, 'ro-', linewidth=2)
                self.ax5.set_xlabel('Principal Component')
                self.ax5.set_ylabel('Proportion of Variance Explained')
                self.canvas4.draw()
         

                ### Make 3d matrix of scores , each plane one PC
                scores_cube = np.zeros((2 , self.size_x , self.size_y))
                i = 0
                for counter5 in range(0,self.size_x):
                    for counter6 in range(0,self.size_y):
                        scores_cube[: , counter5 , counter6] = scores[i, 0:2]
                        i += 1


                
                ##### PC1 and PC2 score value images

                if self.radioButton_pcacolormap.isChecked() == False:
                    pca_cmap = 'jet'
                else:
                    pca_cmap = 'viridis'

                im0 = self.ax6.imshow(scores_cube[0 , : , :], cmap = pca_cmap, interpolation ='nearest')
                divider = make_axes_locatable(self.ax6)
                cax6 = divider.append_axes("right", size="5%", pad=0.05) 
                self.cb1 = self.fig5.colorbar(im0 ,  cax=cax6)
                self.ax6.set_yticks([])
                self.ax6.set_xticks([])
                self.ax6.set_yticklabels([])
                self.ax6.set_xticklabels([])


                im1 = self.ax7.imshow(scores_cube[1 , : , :], cmap = pca_cmap, interpolation ='nearest') 
                divider = make_axes_locatable(self.ax7)
                cax7 = divider.append_axes("right", size="5%", pad=0.05)
                self.cb2 = self.fig5.colorbar(im1 , cax=cax7)
                self.ax7.set_yticks([])
                self.ax7.set_xticks([])
                self.ax7.set_yticklabels([])
                self.ax7.set_xticklabels([])


                self.ax6.set_title('PC 1')
                self.ax7.set_title('PC 2')
                plt.tight_layout()
                self.canvas5.draw()
                self.check_colorbar = 1

                self.label_note_preprocessing_pca.setText(
                     f"Preprocessing: {self.note_unit_conversion} | {self.note_spike_removal} | {self.note_smoothing} | "
                     f"{self.note_baseline} | {self.note_normalization}")
                
                # print(scores[:,0:2])
                # print(scores_cube)
                self.check_pca = 1

    def OpenFileCls(self):

        file_path_open= QFileDialog.getOpenFileName(self, 'Select the file' , filter = "Text files (*.txt *.asc)")
        sender = self.sender()

        if file_path_open[0] == '': # For cancelling 
            pass
        else:
            df= pd.read_table(file_path_open[0] , header=0 , decimal="," , index_col=False , skiprows=34)
            self.x_ref_cls = df.iloc[:,0]
            self.x_ref_cls = self.x_ref_cls.values
            self.x_ref_cls = self.x_ref_cls.astype(np.float)
            y_ref_cls = df.iloc[:,1]
            y_ref_cls = y_ref_cls.values
            y_ref_cls = y_ref_cls.astype(np.float)                        

            ###Preprocessing###
            # y_ref_cls = self.correct_zero (y_ref_cls)
            # y_ref_cls = rp.normalise(y_ref_cls,self.x_ref_cls,method="area")

            # baseline_roi1 = []
            # for i in range (0,4500,500):
            #     a = [i,i+500]
            #     baseline_roi1.append(a)
            # baseline_roi = np.array(baseline_roi1)
            # y_ref_cls, base_arps2 = rp.baseline(self.x_ref_cls , y_ref_cls, \
            # baseline_roi,'arPLS',lam=10**5,ratio=0.01)
            ######

            if sender.text() == 'Reference 1':
                self.check_clsref1 = 1
                self.y_ref_cls_1 = y_ref_cls
                return (self.y_ref_cls_1)
            elif sender.text() == 'Reference 2':
                self.check_clsref2 = 1
                self.y_ref_cls_2 = y_ref_cls
                return (self.y_ref_cls_2)
            elif sender.text() == 'Mixture':
                self.check_clssam = 1
                self.y_sam_cls = y_ref_cls
                return (self.y_sam_cls)

    def perform_cls (self):

        if self.check_clsref1 == 1 & self.check_clsref2 == 1 & self.check_clssam == 1:

            comparision1 = self.y_ref_cls_1 == self.y_ref_cls_2 
            comparision2 = self.y_ref_cls_1 == self.y_sam_cls 
            comparision3 = self.y_ref_cls_1 == self.y_sam_cls

            if comparision1.all() == True or comparision2.all() == True or comparision3.all() == True:

                msg = QMessageBox()
                msg.setWindowTitle("Error")
                msg.setText("Please select the right references or mixture spectra!")
                msg.setIcon(QMessageBox.Critical) 
                msg.setStandardButtons(QMessageBox.Ok)
                x = msg.exec_()
                if x == QMessageBox.Ok:
                    pass

            else:
                self.ax8.clear()


                k = np.transpose(np.column_stack ((self.y_ref_cls_1 , self.y_ref_cls_2)))
                CLS_concentration = np.dot(inv(np.dot(k , k.T)) , np.dot(k , self.y_sam_cls))
                CLS_results = np.array([['Component 1' , float(CLS_concentration[0])], ['Component 2', \
                                                                           float(CLS_concentration[1])]])


                ratio1 = CLS_concentration[0] * 100 / (CLS_concentration[0] + CLS_concentration[1])
                ratio2 = CLS_concentration[1] * 100 / (CLS_concentration[0] + CLS_concentration[1]) 
                
                self.lcdNumber_component1.display(float(ratio1))
                self.lcdNumber_component2.display(float(ratio2))

                self.ax8.plot(self.x_ref_cls , (CLS_concentration[0] * self.y_ref_cls_1))
                self.ax8.plot(self.x_ref_cls , (CLS_concentration[1] * self.y_ref_cls_2))
                self.ax8.plot(self.x_ref_cls , (self.y_sam_cls))
                self.ax8.set_xlabel('Wavenumber ($cm^{-1}$)')
                self.ax8.set_ylabel('Intensity (a.u.)') 
                self.canvas6.draw()
        else:
            pass
        

    def KmeansClustering (self):
        
        if self.check_run_preprocessing == 0:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information) 
            msg.setWindowTitle("Notice")
            msg.setText("Please do preprocessing (preferably with Area normalization), " 
                "Before performing K-Means Cluster Analysis (KMCA).")
            msg.setStandardButtons(QMessageBox.Ok)
            x = msg.exec_()
            if x == QMessageBox.Ok:
                pass
            
        else:
            if self.size_x == 1 and self.size_y == 1:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Information) 
                msg.setWindowTitle("Note")
                msg.setText("For performing K-Means Cluster Analysis, more than one Raman data is needed.")
                msg.setStandardButtons(QMessageBox.Ok)
                x = msg.exec_()
                if x == QMessageBox.Ok:
                    pass

            else: 
                if self.lineEdit_7_kmeans_clusters.text() == '':
                    msg1 = QMessageBox()
                    msg1.setWindowTitle("Error")
                    msg1.setText("Please enter the number of clusters !")
                    msg1.setIcon(QMessageBox.Critical) 
                    msg1_b = msg1.setStandardButtons(QMessageBox.Ok)
                    return_value = msg1.exec_()
                    if return_value == QMessageBox.Ok:
                        pass

                elif int(self.lineEdit_7_kmeans_clusters.text()) <= 1:
                    msg1 = QMessageBox()
                    msg1.setWindowTitle("Error")
                    msg1.setText("The number of clusters must be more than 1 !")
                    msg1.setIcon(QMessageBox.Critical) 
                    msg1_b = msg1.setStandardButtons(QMessageBox.Ok)
                    return_value = msg1.exec_()
                    if return_value == QMessageBox.Ok:
                        pass

                elif int(self.lineEdit_7_kmeans_clusters.text()) > len(self.pdm):
                    msg1 = QMessageBox()
                    msg1.setWindowTitle("Error")
                    msg1.setText("The number of clusters must be equal or less than the number of input Data !")
                    msg1.setIcon(QMessageBox.Critical) 
                    msg1_b = msg1.setStandardButtons(QMessageBox.Ok)
                    return_value = msg1.exec_()
                    if return_value == QMessageBox.Ok:
                        pass
                else:
                    
                    self.num_kmeans_clusters = self.lineEdit_7_kmeans_clusters.text()

                    if self.check_kmeans_clear == 1: # After first run clears the plots
                        self.fig7.clear()
                        self.layout_9.removeWidget(self.toolbar7)
                        self.layout_9.removeWidget(self.canvas7)
                        # self.cb3.remove()

                        self.fig8.clear()
                        self.layout_10.removeWidget(self.toolbar8)
                        self.layout_10.removeWidget(self.canvas8)

                        self.fig9.clear()
                        self.layout_11.removeWidget(self.toolbar9)
                        self.layout_11.removeWidget(self.canvas9)
     
                    self.pdm_kmca = self.pdm  ## For MCR
                    self.x_ref1_new_kmca = self.x_ref1_new ## For MCR

                    self.check_kmeans_clear = 1 

                    # Figure for Image
                    self.fig7 = plt.figure('K-Means Image')
                    self.ax9 = self.fig7.add_subplot(111)
                    self.canvas7 = FigureCanvas(self.fig7)
                    self.toolbar7 = NavigationToolbar(self.canvas7, self)
                    self.layout_9.addWidget(self.toolbar7)
                    self.layout_9.addWidget(self.canvas7)
                    
                    # Figure for plots
                    self.fig8 , self.axs = plt.subplots(int(self.num_kmeans_clusters) , sharex=True)
                    self.canvas8 = FigureCanvas(self.fig8)
                    self.toolbar8 = NavigationToolbar(self.canvas8, self)
                    self.layout_10.addWidget(self.toolbar8)
                    self.layout_10.addWidget(self.canvas8)

                    # Figure for elbow plot
                    self.fig9 = plt.figure('Elbow Plot')
                    self.ax10 = self.fig9.add_subplot(111)
                    self.canvas9 = FigureCanvas(self.fig9)
                    self.toolbar9 = NavigationToolbar(self.canvas9, self)
                    self.layout_11.addWidget(self.toolbar9)
                    self.layout_11.addWidget(self.canvas9)

                    distortions = []
                    K = range(1,int(self.num_kmeans_clusters) + 1)
                    for k in K:                    
                        kmeans = KMeans(n_clusters = k)
                        kmeans.fit(self.pdm)
                        distortions.append(kmeans.inertia_)

                    self.kmeans_cluster_centers = kmeans.cluster_centers_ ### Made for open separate windows
                    self.kmeans_labels = kmeans.labels_

                    ### Make K-Means reconstructed images
                    i = 0
                    self.kmean_matrix = np.zeros((self.size_x,self.size_y))
                    for counter9 in range(0,self.size_x):
                        for counter10 in range(0,self.size_y):
                            self.kmean_matrix[counter9,counter10] = kmeans.labels_[i] + 1
                            i += 1

                    print(self.kmean_matrix)
                    # plt.figure('Kmeans')
                    cmap = plt.cm.get_cmap("gist_rainbow", self.num_kmeans_clusters)
                    im3 = self.ax9.imshow(self.kmean_matrix[:,:], cmap = cmap, interpolation ='nearest')

                    # Modify the colorbar
                    # divider = make_axes_locatable(self.ax9)
                    # cax9 = divider.append_axes("right", size="5%", pad=0.05) 
                    self.cb3 = self.fig7.colorbar(im3 , ax=self.ax9 , \
                        ticks = np.linspace(1,int(self.num_kmeans_clusters), int(self.num_kmeans_clusters)))
                    self.ax9.set_yticks([])
                    self.ax9.set_xticks([])
                    self.ax9.set_yticklabels([])
                    self.ax9.set_xticklabels([])
                    self.canvas7.draw()

                    ### Plot K-Means mean spectra , Save
                    j = 0
                    for j in range(0, len(self.pdm)): # size of PDM
                        self.axs[kmeans.labels_[j]].clear()   
                        self.axs[kmeans.labels_[j]].plot(self.x_ref1_new, kmeans.cluster_centers_[kmeans.labels_[j],:])
                        self.axs[kmeans.labels_[j]].set_title(f'{kmeans.labels_[j] + 1}', y=1.0, pad=-15, fontsize=11, loc='right')
                        if self.check_save_kmeans == 1:
                            df = pd.DataFrame(np.column_stack((self.x_ref1_new, kmeans.cluster_centers_[kmeans.labels_[j],:])))
                            df.to_csv(self.directory_path_save_kmean + '/' + 'New_' + f'K-Means - {kmeans.labels_[j] + 1}' + '.txt',\
                                         sep='\t', header=False, index=False)
                        # print(kmeans.labels_[j])

                    self.canvas8.draw()


                    ### Plot Elbow
                    self.ax10.plot(K, distortions, 'bx-')
                    plt.xlabel('Number of Clusters')
                    plt.ylabel('Distortion')
                    plt.title('The Elbow Method Showing the Optimal Number of Clusters')
                    #Sum of Squared Distances From Each Point to Its Assigned Center
                    self.canvas9.draw()


                    ### Initial Spectra For MCR
                    self.initial_spectra_kmca = np.vstack((kmeans.cluster_centers_[0], kmeans.cluster_centers_[1]))
                    if kmeans.cluster_centers_.shape[0] > 2:
                        for i in range(2,kmeans.cluster_centers_.shape[0]):
                            self.initial_spectra_kmca = np.vstack((self.initial_spectra_kmca, kmeans.cluster_centers_[i]))

                    if self.norm_area_check == 1 :
                        self.norm_area_check_kmeans = 1 
                    elif self.norm_area_check == 0 :
                        self.norm_area_check_kmeans = 0
                    self.check_kmca = 1


                    ### Button for separate windows

                    self.Button_openseparate_kmca = QtWidgets.QPushButton(self.frame_11)
                    font = QtGui.QFont()
                    font.setFamily("Roboto")
                    self.Button_openseparate_kmca.setFont(font)
                    self.Button_openseparate_kmca.setFocusPolicy(QtCore.Qt.NoFocus)
                    self.Button_openseparate_kmca.setObjectName("Button_openseparate_kmca")
                    self.Button_openseparate_kmca.setGeometry(QtCore.QRect(300, 11, 140, 28))
                    self.Button_openseparate_kmca.setText("Open in Separate Windows")
                    self.Button_openseparate_kmca.clicked.connect(self.OpenPlots_kmca)
                    self.Button_openseparate_kmca.show()


                    self.label_note_preprocessing_km.setText(
                     f"Preprocessing: {self.note_unit_conversion} | {self.note_spike_removal} | {self.note_smoothing} | "
                     f"{self.note_baseline} | {self.note_normalization}")


    def perform_mcr (self):

        if self.radioButton_initial_kmeans.isChecked() == True and self.check_kmca == 0:
        # if self.check_kmca == 0 and self.check_nmf == 0 or : #or self.norm_area_check == 0
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information) 
            msg.setWindowTitle("Notice")
            msg.setText("Before performing Multivariate Curve Resolution Analysis (MCR)"
                ", preprocessing, and K-Means Cluster Analysis (KMCA)"
                " or Non-negative Matrix Factorization (NMF) must be done.")
            msg.setStandardButtons(QMessageBox.Ok)
            x = msg.exec_()
            if x == QMessageBox.Ok:
                pass
        
        elif self.radioButton_initial_nmf.isChecked() == True and self.check_nmf == 0:    
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information) 
            msg.setWindowTitle("Notice")
            msg.setText("Before performing Multivariate Curve Resolution Analysis (MCR)"
                ", preprocessing, and K-Means Cluster Analysis (KMCA)"
                " or Non-negative Matrix Factorization (NMF) must be done.")
            msg.setStandardButtons(QMessageBox.Ok)
            x = msg.exec_()
            if x == QMessageBox.Ok:
                pass

        else:

            if self.radioButton_initial_kmeans.isChecked() == True:
                if self.norm_area_check_kmeans == 0 :
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Information) 
                    msg.setWindowTitle("Notice")
                    msg.setText("Preprocessing with Area normalization (before K-Means cluster analysis) yields better MCR-ALS results.")
                    msg.setInformativeText("Do: Proprocessing with Area normalization -> KMCA -> MCR-ALS")
                    msg.setStandardButtons(QMessageBox.Ok)
                    x = msg.exec_()
                    if x == QMessageBox.Ok:
                        pass

            if self.check_mcr_clear == 1:

                self.fig10.clear()
                self.layout_13.removeWidget(self.toolbar10)
                self.layout_13.removeWidget(self.canvas10)
                # self.cb3.remove()

                self.fig11.clear()
                self.layout_12.removeWidget(self.toolbar11)
                self.layout_12.removeWidget(self.canvas11)


            self.check_mcr_clear = 1 

            ###### MCR-NNLS
            
            if self.radioButton_initial_kmeans.isChecked() == True:
                self.initial_spectra = self.initial_spectra_kmca
                mcrar = McrAR(max_iter=1000, st_regr='NNLS', c_regr='NNLS',\
                    st_constraints=[ConstraintNonneg()],
                    c_constraints=[ConstraintNonneg(), ConstraintNorm()])
                mcrar.fit(self.pdm_kmca, ST= self.initial_spectra)
                self.x_mcr = self.x_ref1_new_kmca 
            
            elif self.radioButton_initial_nmf.isChecked() == True:
                self.initial_spectra = self.nmf_comp
                mcrar = McrAR(max_iter=1000, st_regr='NNLS', c_regr='NNLS',\
                    st_constraints=[ConstraintNonneg()],
                    c_constraints=[ConstraintNonneg(), ConstraintNorm()]) #, ConstraintNorm()
                # con_ = np.array([[0 , 1] , [0.2 , 0.8] , [0.4 , 0.6] , [1 , 0] , [0.8 , 0.2] , [0.6 , 0.4]])
                mcrar.fit(self.pdm_nmf, ST= self.initial_spectra) 
                # C = con_ , c_fix=[0,1], st_fix=[0] , c_first=True
                self.x_mcr = self.x_ref1_new_nmf 

            # Make Figures for Images
            if int(self.initial_spectra.shape[0]) <= 3:
                self.fig10 , axs2 = plt.subplots(1 , self.initial_spectra.shape[0] , sharex=True)
                self.fig10.tight_layout(pad=1.25)

            elif int(self.initial_spectra.shape[0]) > 3 and int(self.initial_spectra.shape[0]) <= 6:
                self.fig10 , axs2 = plt.subplots(2 , math.ceil(self.initial_spectra.shape[0]/2) , sharex=True)
                self.fig10.tight_layout(pad=1.25)
                if self.initial_spectra.shape[0] % 2 != 0:
                    axs2[1,-1].patch.set_visible(False)
                    axs2[1,-1].axis('off')

            elif int(self.initial_spectra.shape[0]) > 6:
                self.fig10 , axs2 = plt.subplots(3 , math.ceil(self.initial_spectra.shape[0]/3) , sharex=True)
                self.fig10.tight_layout(pad=1.25)
                if self.initial_spectra.shape[0] % 3 == 1:
                    axs2[2,-1].patch.set_visible(False)
                    axs2[2,-1].axis('off')
                    axs2[2,-2].patch.set_visible(False)
                    axs2[2,-2].axis('off')
                elif self.initial_spectra.shape[0] % 3 == 2:
                    axs2[2,-1].patch.set_visible(False)
                    axs2[2,-1].axis('off')
                    

            self.canvas10 = FigureCanvas(self.fig10)
            self.toolbar10 = NavigationToolbar(self.canvas10, self)
            self.layout_13.addWidget(self.toolbar10)
            self.layout_13.addWidget(self.canvas10)
            
            # Figure for plots
            self.fig11 , axs3 = plt.subplots(self.initial_spectra.shape[0] , sharex=True)
            self.canvas11 = FigureCanvas(self.fig11)
            self.toolbar11 = NavigationToolbar(self.canvas11, self)
            self.layout_12.addWidget(self.toolbar11)
            self.layout_12.addWidget(self.canvas11)

            
            self.mcrar_shape = mcrar.C_opt_.shape
            self.mcrar_ST_opt = mcrar.ST_opt_

            #### Plot Spectra , Save
            for j in range(mcrar.C_opt_.shape[1]):
                axs3[j].plot(self.x_mcr, mcrar.ST_opt_[j,:])
                axs3[j].set_title(f'{j+1}', y=1.0, pad = -15, fontsize=11, loc='right')
                if self.check_save_mcr == 1:
                    df = pd.DataFrame(np.column_stack((self.x_mcr, mcrar.ST_opt_[j,:])))
                    df.to_csv(self.directory_path_save_mcr + '/' + 'New_' + f'MCR-ALS - {j+1}' + '.txt',\
                                         sep='\t', header=False, index=False)


            self.canvas11.draw()
            self.Button_openseparate_mcr = QtWidgets.QPushButton(self.frame_14)
            font = QtGui.QFont()
            font.setFamily("Roboto")
            self.Button_openseparate_mcr.setFont(font)
            self.Button_openseparate_mcr.setFocusPolicy(QtCore.Qt.NoFocus)
            self.Button_openseparate_mcr.setObjectName("Button_openseparate_kmca")
            self.Button_openseparate_mcr.setGeometry(QtCore.QRect(300, 11, 140, 28))
            self.Button_openseparate_mcr.setText("Open in Separate Windows")
            self.Button_openseparate_mcr.clicked.connect(self.OpenPlots_mcr)
            self.Button_openseparate_mcr.show()


            # mcrar.C_opt_ = minmax_scale(mcrar.C_opt_, feature_range=(0,1)) # To change the scale to 0 and 1
            
            # mcrar.C_opt_[:,0] = minmax_scale(mcrar.C_opt_[:,0], feature_range=(0,1)) # To change the range of each component to a certain number
            # mcrar.C_opt_[:,1] = minmax_scale(mcrar.C_opt_[:,1], feature_range=(0,1))
            # mcrar.C_opt_[:,2] = minmax_scale(mcrar.C_opt_[:,2], feature_range=(0,1))

            #### MCR Images
            t = 0
            mcr_matrix = np.zeros((mcrar.C_opt_.shape[1], self.size_x,self.size_y))
            for d in  range(0 , mcrar.C_opt_.shape[1]):
                for a in range(0, self.size_x):
                    for b in range(0, self.size_y):
                        mcr_matrix[d, a , b] = mcrar.C_opt_[t,d]
                        t += 1
                t = 0

            t = 0
            for j in range(mcrar.C_opt_.shape[1]): ### Make the images and put them on a location in screen, considering the number of images (components)
                if mcrar.C_opt_.shape[1] <= 3:
                    im = axs2[j].imshow(mcr_matrix[j , : , :] * 100, cmap = 'viridis', interpolation ='nearest',vmin=0, vmax=100)
                    axs2[j].set_title(f'{j+1}', fontsize=11, loc='left')
                    axs2[j].set_yticks([])
                    axs2[j].set_xticks([])
                    axs2[j].set_yticklabels([])
                    axs2[j].set_xticklabels([])
                    divider = make_axes_locatable(axs2[j])
                    cax = divider.append_axes("right", size="5%", pad=0.05)
                    self.fig10.colorbar(im , cax = cax)
                    cbar = self.fig10.colorbar(im , cax = cax)

                    

                elif mcrar.C_opt_.shape[1] > 3 and mcrar.C_opt_.shape[1] <= 6:
                    if j + 1 <= math.ceil(mcrar.C_opt_.shape[1]/2):    
                        im = axs2[0,j].imshow(mcr_matrix[j , : , :], cmap = 'viridis', interpolation ='nearest',vmin=0, vmax=100)
                        axs2[0,j].set_title(f'{j+1}', fontsize=11, loc='left')
                        axs2[0,j].set_yticks([])
                        axs2[0,j].set_xticks([])
                        axs2[0,j].set_yticklabels([])
                        axs2[0,j].set_xticklabels([])
                        divider = make_axes_locatable(axs2[0,j])
                        cax = divider.append_axes("right", size="5%", pad=0.05)
                        cbar = self.fig10.colorbar(im , cax = cax)
                        
                    if j + 1 > math.ceil(mcrar.C_opt_.shape[1]/2):    
                        im = axs2[1,j - math.ceil(mcrar.C_opt_.shape[1]/2)].imshow(mcr_matrix[j , : , :], cmap = 'viridis', interpolation ='nearest')
                        axs2[1,j - math.ceil(mcrar.C_opt_.shape[1]/2)].set_title(f'{j+1}', fontsize=11, loc='left')
                        axs2[1,j - math.ceil(mcrar.C_opt_.shape[1]/2)].set_yticks([])
                        axs2[1,j - math.ceil(mcrar.C_opt_.shape[1]/2)].set_xticks([])
                        axs2[1,j - math.ceil(mcrar.C_opt_.shape[1]/2)].set_yticklabels([])
                        axs2[1,j - math.ceil(mcrar.C_opt_.shape[1]/2)].set_xticklabels([])
                        divider = make_axes_locatable(axs2[1,j - math.ceil(mcrar.C_opt_.shape[1]/2)])
                        cax = divider.append_axes("right", size="5%", pad=0.05)
                        self.fig10.colorbar(im , cax = cax)
                        cbar = self.fig10.colorbar(im , cax = cax)
                         
                elif mcrar.C_opt_.shape[1] >= 7:
                    if j + 1 <= math.ceil(mcrar.C_opt_.shape[1]/3):    
                        im = axs2[0,j].imshow(mcr_matrix[j , : , :], cmap = 'viridis', interpolation ='nearest',vmin=0, vmax=100)
                        axs2[0,j].set_title(f'{j+1}', fontsize=11, loc='left')
                        axs2[0,j].set_yticks([])
                        axs2[0,j].set_xticks([])
                        axs2[0,j].set_yticklabels([])
                        axs2[0,j].set_xticklabels([])
                        divider = make_axes_locatable(axs2[0,j])
                        cax = divider.append_axes("right", size="5%", pad=0.05)
                        cbar = self.fig10.colorbar(im , cax = cax)
                    
                    if j + 1 > math.ceil(mcrar.C_opt_.shape[1]/3) and j + 1 <= math.ceil(mcrar.C_opt_.shape[1]/3) * 2:    
                        im = axs2[1,j - math.ceil(mcrar.C_opt_.shape[1]/3)].imshow(mcr_matrix[j , : , :], cmap = 'viridis', interpolation ='nearest')
                        axs2[1,j - math.ceil(mcrar.C_opt_.shape[1]/3)].set_title(f'{j+1}', fontsize=11, loc='left')
                        axs2[1,j - math.ceil(mcrar.C_opt_.shape[1]/3)].set_yticks([])
                        axs2[1,j - math.ceil(mcrar.C_opt_.shape[1]/3)].set_xticks([])
                        axs2[1,j - math.ceil(mcrar.C_opt_.shape[1]/3)].set_yticklabels([])
                        axs2[1,j - math.ceil(mcrar.C_opt_.shape[1]/3)].set_xticklabels([])
                        divider = make_axes_locatable(axs2[1,j - math.ceil(mcrar.C_opt_.shape[1]/3)])
                        cax = divider.append_axes("right", size="5%", pad=0.05)
                        self.fig10.colorbar(im , cax = cax)
                        cbar = self.fig10.colorbar(im , cax = cax)

                    if j + 1 > math.ceil((mcrar.C_opt_.shape[1]/3)) * 2:    
                        im = axs2[2,j - math.ceil(mcrar.C_opt_.shape[1]/3) * 2].imshow(mcr_matrix[j , : , :], cmap = 'viridis', interpolation ='nearest')
                        axs2[2,j - math.ceil(mcrar.C_opt_.shape[1]/3) * 2].set_title(f'{j+1}', fontsize=11, loc='left')
                        axs2[2,j - math.ceil(mcrar.C_opt_.shape[1]/3) * 2].set_yticks([])
                        axs2[2,j - math.ceil(mcrar.C_opt_.shape[1]/3) * 2].set_xticks([])
                        axs2[2,j - math.ceil(mcrar.C_opt_.shape[1]/3) * 2].set_yticklabels([])
                        axs2[2,j - math.ceil(mcrar.C_opt_.shape[1]/3) * 2].set_xticklabels([])
                        divider = make_axes_locatable(axs2[2,j - math.ceil(mcrar.C_opt_.shape[1]/3) * 2])
                        cax = divider.append_axes("right", size="5%", pad=0.05)
                        self.fig10.colorbar(im , cax = cax)
                        cbar = self.fig10.colorbar(im , cax = cax)

            print(mcrar.C_opt_)
            self.canvas10.draw()

            
    def perform_nmf(self):
        if self.check_run_preprocessing == 0:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information) 
            msg.setWindowTitle("Notice")
            msg.setText("Please do preprocessing before Performing Non-negative Matrix Factorization (NMF).")
            msg.setStandardButtons(QMessageBox.Ok)
            x = msg.exec_()
            if x == QMessageBox.Ok:
                pass
        else:
            if self.size_x == 1 and self.size_y == 1:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Information) 
                msg.setWindowTitle("Note")
                msg.setText("For performing Non-negative Matrix Factorization, more than one Raman data is needed.")
                msg.setStandardButtons(QMessageBox.Ok)
                x = msg.exec_()
                if x == QMessageBox.Ok:
                    pass 

            else:                                
                self.num_nmf_comp = self.lineEdit_8_nmf_num.text()

                if self.num_nmf_comp == '' or int(self.num_nmf_comp) < 2:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Critical) 
                    msg.setWindowTitle("Notice")
                    if self.num_nmf_comp == '':
                        msg.setText("Please enter the number of components!")
                    elif int(self.num_nmf_comp) < 2:
                        msg.setText("The number of components must be larger than 1 !")
                    msg.setStandardButtons(QMessageBox.Ok)
                    x = msg.exec_()
                    if x == QMessageBox.Ok:
                        pass
                else:
                    if self.check_nmf_clear == 1: # After first run clears the plots
                        self.fig12.clear()
                        self.layout_14.removeWidget(self.toolbar12)
                        self.layout_14.removeWidget(self.canvas12)

                    self.pdm_nmf = self.pdm   ## For MCR
                    self.x_ref1_new_nmf = self.x_ref1_new ## For MCR

                    self.check_nmf_clear = 1
                    self.fig12 , axs4 = plt.subplots(int(self.num_nmf_comp) , sharex=True)
                    self.canvas12 = FigureCanvas(self.fig12)
                    self.toolbar12 = NavigationToolbar(self.canvas12, self)
                    self.layout_14.addWidget(self.toolbar12)
                    self.layout_14.addWidget(self.canvas12)


                    self.pdm[self.pdm < 0] = 0
                    nmf = decomposition.NMF(n_components = int(self.num_nmf_comp), random_state=None , max_iter=500)
                    
                    nmf.fit(self.pdm)
                    self.nmf_comp = nmf.components_
                    j = 0
                    for j in range(0, int(self.num_nmf_comp)): # number of components
                        axs4[j].clear()   
                        axs4[j].plot(self.x_ref1_new, self.nmf_comp[j,:])
                        axs4[j].set_title(f'{j + 1}', y=1.0, pad=-15, fontsize=11, loc='right')
                        if self.check_save_nmf == 1:
                            df = pd.DataFrame(np.column_stack((self.x_ref1_new, self.nmf_comp[j,:])))
                            df.to_csv(self.directory_path_save_nmf + '/' + 'New_' + f'NMF - {j + 1}' + '.txt',\
                                         sep='\t', header=False, index=False)
                        
                    self.check_nmf = 1
                    self.canvas12.draw()
                
                    self.Button_openseparate_nmf = QtWidgets.QPushButton(self.frame_17)
                    font = QtGui.QFont()
                    font.setFamily("Roboto")
                    self.Button_openseparate_nmf.setFont(font)
                    self.Button_openseparate_nmf.setFocusPolicy(QtCore.Qt.NoFocus)
                    self.Button_openseparate_nmf.setObjectName("Button_openseparate_nmf")
                    self.Button_openseparate_nmf.setGeometry(QtCore.QRect(300, 11, 140, 28))
                    self.Button_openseparate_nmf.setText("Open in Separate Windows")
                    self.Button_openseparate_nmf.clicked.connect(self.OpenPlots_nmf)
                    self.Button_openseparate_nmf.show()


                    self.label_note_preprocessing_nmf.setText(
                     f"Preprocessing: {self.note_unit_conversion} | {self.note_spike_removal} | {self.note_smoothing} | "
                     f"{self.note_baseline} | {self.note_normalization}")


            
            

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
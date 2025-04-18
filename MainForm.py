# Form implementation generated from reading ui file 'form.ui'
#
# Created by: PyQt6 UI code generator 6.8.1
#
# WARNING: Any manual changes made to this file will be lost when pyuic6 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QMessageBox  # Import QMessageBox
from draw import Draw
from algorithms import *

class Ui_MainForm(object):
    def setupUi(self, MainForm):
        # set up the main window
        MainForm.setObjectName("MainForm")
        MainForm.resize(1755, 1327)
        self.centralwidget = QtWidgets.QWidget(parent=MainForm)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.Canvas = Draw(parent=self.centralwidget)
        self.Canvas.setObjectName("Canvas")
        self.horizontalLayout.addWidget(self.Canvas)
        MainForm.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(parent=MainForm)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1755, 22))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(parent=self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuSimplify = QtWidgets.QMenu(parent=self.menubar)
        self.menuSimplify.setObjectName("menuSimplify")
        self.menuView = QtWidgets.QMenu(parent=self.menubar)
        self.menuView.setObjectName("menuView")
        MainForm.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(parent=MainForm)
        self.statusbar.setObjectName("statusbar")
        MainForm.setStatusBar(self.statusbar)
        self.toolBar = QtWidgets.QToolBar(parent=MainForm)
        self.toolBar.setObjectName("toolBar")
        MainForm.addToolBar(QtCore.Qt.ToolBarArea.TopToolBarArea, self.toolBar)
        self.actionOpen = QtGui.QAction(parent=MainForm)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("images/icons/open_file.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionOpen.setIcon(icon)
        self.actionOpen.setObjectName("actionOpen")
        self.actionExit = QtGui.QAction(parent=MainForm)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap("images/icons/exit.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionExit.setIcon(icon1)
        self.actionExit.setObjectName("actionExit")
        self.actionMBR = QtGui.QAction(parent=MainForm)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap("images/icons/maer.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionMBR.setIcon(icon2)
        self.actionMBR.setObjectName("actionMBR")
        self.actionPCA = QtGui.QAction(parent=MainForm)
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap("images/icons/pca.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionPCA.setIcon(icon3)
        self.actionPCA.setObjectName("actionPCA")
        self.actionClear_results = QtGui.QAction(parent=MainForm)
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap("images/icons/clear.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionClear_results.setIcon(icon4)
        self.actionClear_results.setObjectName("actionClear_results")
        self.actionClear_all = QtGui.QAction(parent=MainForm)
        icon5 = QtGui.QIcon()
        icon5.addPixmap(QtGui.QPixmap("images/icons/clear_er.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionClear_all.setIcon(icon5)
        self.actionClear_all.setObjectName("actionClear_all")
        self.actionLongestEdge = QtGui.QAction(parent=MainForm)
        icon6 = QtGui.QIcon()
        icon6.addPixmap(QtGui.QPixmap("images/icons/longestedge.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionLongestEdge.setIcon(icon6)
        self.actionLongestEdge.setObjectName("actionLongestEdge")
        self.actionWallAverage = QtGui.QAction(parent=MainForm)
        icon7 = QtGui.QIcon()
        icon7.addPixmap(QtGui.QPixmap("images/icons/wa.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionWallAverage.setIcon(icon7)
        self.actionWallAverage.setObjectName("actionWallAverage")
        self.actionWeightedBisector = QtGui.QAction(parent=MainForm)
        icon8 = QtGui.QIcon()
        icon8.addPixmap(QtGui.QPixmap("images/icons/weightedbisector.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionWeightedBisector.setIcon(icon8)
        self.actionWeightedBisector.setObjectName("actionWeightedBisector")
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionExit)
        self.menuSimplify.addAction(self.actionMBR)
        self.menuSimplify.addAction(self.actionPCA)
        self.menuSimplify.addAction(self.actionLongestEdge)
        self.menuSimplify.addAction(self.actionWallAverage)
        self.menuSimplify.addAction(self.actionWeightedBisector)
        self.menuView.addAction(self.actionClear_results)
        self.menuView.addAction(self.actionClear_all)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuSimplify.menuAction())
        self.menubar.addAction(self.menuView.menuAction())
        self.toolBar.addAction(self.actionOpen)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionMBR)
        self.toolBar.addAction(self.actionPCA)
        self.toolBar.addAction(self.actionLongestEdge)
        self.toolBar.addAction(self.actionWallAverage)
        self.toolBar.addAction(self.actionWeightedBisector)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionClear_results)
        self.toolBar.addAction(self.actionClear_all)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionExit)

        self.retranslateUi(MainForm)
        QtCore.QMetaObject.connectSlotsByName(MainForm)
        
        #Connect menu and functions
        self.actionOpen.triggered.connect(self.openShapefile)
        self.actionMBR.triggered.connect(self.simplifyBuildingMBR)
        self.actionPCA.triggered.connect(self.simplifyBuildingBRPCA)
        self.actionLongestEdge.triggered.connect(self.simplifyBuildingLongestEdge)
        self.actionWallAverage.triggered.connect(self.simplifyBuildingWallAverage)
        self.actionWeightedBisector.triggered.connect(self.simplifyBuildingWeightedBisector)
        self.actionClear_results.triggered.connect(self.clearResults)
        self.actionClear_all.triggered.connect(self.clearAll)
        self.actionExit.triggered.connect(self.exitApplication)

    def retranslateUi(self, MainForm):
        '''
        This function sets the text and tooltips for the UI elements.
        
        '''
        _translate = QtCore.QCoreApplication.translate
        MainForm.setWindowTitle(_translate("MainForm", "Building simplify"))
        self.menuFile.setTitle(_translate("MainForm", "File"))
        self.menuSimplify.setTitle(_translate("MainForm", "Simplify"))
        self.menuView.setTitle(_translate("MainForm", "View"))
        self.toolBar.setWindowTitle(_translate("MainForm", "toolBar"))
        self.actionOpen.setText(_translate("MainForm", "Open"))
        self.actionOpen.setToolTip(_translate("MainForm", "Open file"))
        self.actionExit.setText(_translate("MainForm", "Exit"))
        self.actionExit.setToolTip(_translate("MainForm", "Close application"))
        self.actionMBR.setText(_translate("MainForm", "MBR"))
        self.actionMBR.setToolTip(_translate("MainForm", "Simplify building using MBR"))
        self.actionPCA.setText(_translate("MainForm", "PCA"))
        self.actionPCA.setToolTip(_translate("MainForm", "Simplify building using PCA"))
        self.actionClear_results.setText(_translate("MainForm", "Clear results"))
        self.actionClear_all.setText(_translate("MainForm", "Clear all"))
        self.actionClear_all.setToolTip(_translate("MainForm", "Clear all data"))
        self.actionLongestEdge.setText(_translate("MainForm", "Longest Edge"))
        self.actionLongestEdge.setToolTip(_translate("MainForm", "Simplify building using Longest Edge"))
        self.actionWallAverage.setText(_translate("MainForm", "Wall Average"))
        self.actionWallAverage.setToolTip(_translate("MainForm", "Simplify building using Wall Average"))
        self.actionWeightedBisector.setText(_translate("MainForm", "Weighted Bisector"))
        self.actionWeightedBisector.setToolTip(_translate("MainForm", "Simplify building using Weighted Bisector"))


    def openShapefile(self):
        """Open and load shapefile data"""
        result = self.Canvas.loadData()
        if result:
            self.statusbar.showMessage("Shapefile loaded successfully", 5000)
        else:
            self.statusbar.showMessage("Failed to load shapefile", 5000)

    def _show_accuracy_popup(self, method_name, accuracy1, accuracy2, percentage1, percentage2):
        """Helper function to display the accuracy metrics (Δσ₁, Δσ₂, their percentages."""
        msg = QMessageBox()
        msg.setWindowTitle("Generalization Accuracy")

        accuracy1_deg = accuracy1 * 180 / pi
        accuracy2_deg = accuracy2 * 180 / pi
        # Format accuracy values (degrees) and percentages
        msg.setText(f"Method: {method_name}\n"
                    f"Average Δσ₁: {accuracy1_deg:.2f}° ({percentage1:.1f}%)\n"
                    f"Average Δσ₂: {accuracy2_deg:.2f}° ({percentage2:.1f}%)")
        msg.setIcon(QMessageBox.Icon.Information)
        msg.exec()

    def _simplify_building(self, method_name, algorithm_func, generalize_method_key=None):
        """universal building simplification function
        
        Arguments:
            method_name: Display name of the method
            algorithm_func: Function reference to the specific algorithm
            generalize_method_key: Key used in generalize_pols 
        """
        # get input data
        building = self.Canvas.getBuilding()
        
        # gse method_name as the generalize key if not specified
        if generalize_method_key is None:
            generalize_method_key = method_name
            
        if not building.isEmpty():
            # s implify building
            a = Algorithms() # create instance of Algorithms
            building_simp = algorithm_func(a, building) 
            
            # Set results
            self.Canvas.setSimplifBuilding(building_simp)
        
        # simplify loaded polygons if any
        polygons = self.Canvas.getPolygons()
        if polygons:
            a = Algorithms()
            # get simplified polygons + all accuracy metrics
            simplified, avg_accuracy_1, avg_accuracy_2, avg_percentage_1, avg_percentage_2 = a.generalize_pols(polygons, generalize_method_key)
            self.Canvas.setSimplifiedPolygons(simplified)
            # show accuracy popup 
            self._show_accuracy_popup(method_name, avg_accuracy_1, avg_accuracy_2, avg_percentage_1, avg_percentage_2)
        else:
            if not building.isEmpty():
                QMessageBox.information(self.Canvas, "Info", 
                                       "Simplified single building. Accuracy calculation requires loaded polygons.")
            else:
                QMessageBox.warning(self.Canvas, "Warning", 
                                   "No building or polygons loaded to simplify.")
        
    
        self.Canvas.repaint()

    # helper functions that call the universal function with targer methods
    def simplifyBuildingMBR(self):
        self._simplify_building("MBR", lambda a, b: a.createMBR(b))
        
    def simplifyBuildingBRPCA(self):
        self._simplify_building("PCA", lambda a, b: a.createBRPCA(b))
    
    def simplifyBuildingLongestEdge(self):
        self._simplify_building("Longest Edge", lambda a, b: a.longestEdge(b), "LongestEdge")
    
    def simplifyBuildingWallAverage(self):
        self._simplify_building("Wall Average", lambda a, b: a.wallAverage(b), "WallAverage")
    
    def simplifyBuildingWeightedBisector(self):
        self._simplify_building("Weighted Bisector", lambda a, b: a.weighted_bisector(b), "WeightedBisector")

    def clearResults(self):
        """Clear only the simplified results"""
        self.Canvas.clearResults()
        self.statusbar.showMessage("Results cleared", 3000)
    
    def clearAll(self):
        """Clear everything"""
        self.Canvas.clearAllData()
        self.statusbar.showMessage("All data cleared", 3000)
    
    def exitApplication(self):
        """Exit the application"""
        QtWidgets.QApplication.quit()
        

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainForm = QtWidgets.QMainWindow()
    ui = Ui_MainForm()
    ui.setupUi(MainForm)
    MainForm.show()
    sys.exit(app.exec())

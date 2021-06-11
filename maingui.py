from PyQt5 import QtWidgets, uic, Qt
from PyQt5.QtWidgets import QInputDialog, QFileDialog, QTableWidget, QTableWidgetItem, QMessageBox
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import sys, os
import numpy as np
from openbabel import openbabel as ob
atom_radius_lib={'Fe':0.5, 'N':0.3,'P':0.4,'C':0.3,'H':0.15, 'Pt':0.4, 'O':0.25, 'F':0.2, 'S':0.3}
bond_distance_lib={'Fe':2.04, 'N':2, 'C':1.55, 'P':1.9, 'H':1.2, 'Pt':2.2}
atom_type_to_color={'Fe':'orange', 'N':'blue', 'C':'gray', 'P':'yellow', 'H':'lightGray', 'Pt':'white',\
'O':'red', 'F':'lightBlue', 'S':'brown'}


object_map = {}
atom_list = []
atom_type = []
atom_Num = [] 
x = []
y = []
z = []

l_object_map = {} #probably a less wordy/messy way to do this, but these are for the ligand window
l_atom_list = []
l_atom_type = []
l_atom_Num = []
l_x = []
l_y = []
l_z = []


class Atom():  
    def __init__(self, atom_type, x_coord, y_coord, z_coord, radius, m, atom_color):
        self.atom_type = atom_type
        self.x_coord = x_coord
        self.y_coord = y_coord
        self.z_coord = z_coord
        self.radius = radius
        self.m = m
        self.atom_color = atom_color


class MyView(pg.opengl.GLViewWidget):
    def mousePressEvent(self, ev):
        #if selctMultipleButton isn't clicked, then only 1 atom at a time is highlighted
        #if selectionTool == None: 
        #self.unHighLight(None) #we want to unhighlight every time something is clicked on??
        #...although maybe we want the user to click on something again to unhighlight??
        self.mousePos = ev.pos()
        try: #if this try fails, that means the click was in the ligand window
            for i in range(len(x)): #next two lines unhighlight all atoms
                atom_list[i].m.setColor(QtGui.QColor(atom_list[i].atom_color))
            if len(self.itemsAt((ev.pos().x(), ev.pos().y(), 10, 10))) != 0:
                index =  object_map[self.itemsAt((ev.pos().x(), ev.pos().y(), 10, 10))[0]]
                atom_type_edit.setText(atom_list[index].atom_type)
                atomic_num.setText(str(index + 1)) #atom numbering starts at 1, but python starts at 0, so add 1
                x_line_edit.setText(str(atom_list[index].x_coord))
                y_line_edit.setText(str(atom_list[index].y_coord))
                z_line_edit.setText(str(atom_list[index].z_coord))
                #print(temp_atom_type,temp_atomic_num,x1,y1,z1)
                atom_list[index].m.setColor(QtGui.QColor('cyan'))
                ev.accept()  
                ligandOrMain = 'main'
        except KeyError:
            for i in range(len(l_x)): #next two lines unhighlight all atoms
                l_atom_list[i].m.setColor(QtGui.QColor(l_atom_list[i].atom_color))
            if len(self.itemsAt((ev.pos().x(), ev.pos().y(), 10, 10))) != 0:
                index =  l_object_map[self.itemsAt((ev.pos().x(), ev.pos().y(), 10, 10))[0]]
                atom_type_edit.setText(l_atom_list[index].atom_type)
                atomic_num.setText(str(index + 1)) #atom numbering starts at 1, but python starts at 0, so add 1
                x_line_edit.setText(str(l_atom_list[index].x_coord))
                y_line_edit.setText(str(l_atom_list[index].y_coord))
                z_line_edit.setText(str(l_atom_list[index].z_coord))
                #print(temp_atom_type,temp_atomic_num,x1,y1,z1)
                l_atom_list[index].m.setColor(QtGui.QColor('cyan'))
                ev.accept()
                ligandOrMain = 'ligand'
            
    def unHighLight(self,ligandOrMain=''):
        if ligandOrMain == 'main':
            for i in range(len(x)):
                atom_list[i].m.setColor(QtGui.QColor(atom_list[i].atom_color))
        elif ligandOrMain:
            for i in range(len(l_x)):
                l_atom_list[i].m.setColor(QtGui.QColor(l_atom_list[i].atom_color))


def calcDistance(a,b):
    return np.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)
def calcLinesToDraw(self,x ,y ,z, atom_type, atom_Num, mainOrLigand):
    if(mainOrLigand=='main'):
        for i in range(len(x)):
            for q in range(i,len(x)):
                dist = calcDistance((x[i],y[i],z[i]),(x[q],y[q],z[q]))
                #if dist < bond_distance_lib[atom_type[i]] or dist < bond_distance_lib[atom_type[q]]:
                #if dist < 2*ob.GetCovalentRad(atom_Num[i]) or dist < 2*ob.GetCovalentRad(atom_Num[q]):
                if dist < (ob.GetCovalentRad(atom_Num[i]) + ob.GetCovalentRad(atom_Num[q])) * 1.1:
                    x1 = (x[i] + ((x[q]-x[i])* atom_radius_lib[str(atom_type[i])]/2),y[i] + ((y[q]-y[i])* atom_radius_lib[str(atom_type[i])]/2)            ,z[i] + ((z[q]-z[i])* atom_radius_lib[str(atom_type[i])]/2))
                    y1 = (x[q] + ((x[i]-x[q])* atom_radius_lib[str(atom_type[q])]/2),y[q] + ((y[i]-y[q])* atom_radius_lib[str(atom_type[q])]/2),z[q] + ((z[i]-z[q])* atom_radius_lib[str(atom_type[q])]/2))
                    pts = np.array([x1,y1])
                    sh1 = gl.GLLinePlotItem(pos = pts, width = 2.5)
                    sh1.setGLOptions('opaque')
                    self.molWindow.addItem(sh1)  
    else:
        for i in range(len(l_x)):
            for q in range(i,len(l_x)):
                dist = calcDistance((l_x[i],l_y[i],l_z[i]),(l_x[q],l_y[q],l_z[q]))
                #if dist < bond_distance_lib[l_atom_type[i]] or dist < bond_distance_lib[l_atom_type[q]]:
                if dist < 2*ob.GetCovalentRad(l_atom_Num[i]) or dist < 2*ob.GetCovalentRad(l_atom_Num[q]):
                    x1 = (l_x[i] + ((l_x[q]-l_x[i])* atom_radius_lib[str(l_atom_type[i])]/2),l_y[i] + ((l_y[q]-l_y[i])* atom_radius_lib[str(l_atom_type[i])]/2)        ,l_z[i] + ((l_z[q]-l_z[i])* atom_radius_lib[str(l_atom_type[i])]/2))
                    y1 = (l_x[q] + ((l_x[i]-l_x[q])* atom_radius_lib[str(l_atom_type[q])]/2),l_y[q] + ((l_y[i]-l_y[q])* atom_radius_lib[str(l_atom_type[q])]/2),l_z[q] + ((l_z[i]-l_z[q])* atom_radius_lib[str(l_atom_type[q])]/2))
                    pts = np.array([x1,y1])
                    sh1 = gl.GLLinePlotItem(pos = pts, width = 2.5)
                    sh1.setGLOptions('opaque')
                    self.ligWindow.addItem(sh1) #testing ligand window       

#loads the tree structure from the library folder (Libs)
def load_project_structure(self,startpath, tree):
    import os
    from PyQt5.QtWidgets import QTreeWidgetItem
    from PyQt5.QtGui import QIcon
    for element in os.listdir(startpath):
        path_info = startpath + '/' + element
        parent_itm = QTreeWidgetItem(tree, [os.path.basename(element)])
        #this next line allows the path to accessed through Qt.UserRole
        parent_itm.setData(0, QtCore.Qt.UserRole, path_info) 
        if os.path.isdir(path_info):
            parent_itm.setFlags(parent_itm.flags() | QtCore.Qt.ItemIsTristate | QtCore.Qt.ItemIsUserCheckable)
            load_project_structure(self,path_info, parent_itm)
            parent_itm.setIcon(0, QIcon('Assets/folder.png'))
        else:
            parent_itm.setFlags(parent_itm.flags() | QtCore.Qt.ItemIsUserCheckable)
            parent_itm.setIcon(0, QIcon('Assets/file.png'))
            parent_itm.setCheckState(0,QtCore.Qt.Unchecked)
       


class MainWindow(QtWidgets.QMainWindow):
        
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        #Load the UI Page
        uic.loadUi('mainwindow.ui', self)
        #Binding Elements
        self.molWindow = self.findChild(gl.GLViewWidget,'MoleculeWidget')
        self.constraint_form_layout = self.findChild(QtWidgets.QFormLayout,'formLayout_3')
        self.molWindow = MyView(self.molWindow)
        self.molWindow.setGeometry(10,10,800,800)
        self.molWindow.setCameraPosition(distance=20)
        self.projectTreeWidget = self.findChild(QtWidgets.QTreeWidget, 'fileWidget')
        self.projectTreeWidget.itemClicked.connect(self.treeItemClicked)
        self.loadMoleculeButton = self.findChild(QtWidgets.QAction, 'actionLoad_Molecule')
        self.loadMoleculeButton.triggered.connect(self.loadMolecule)
        #ligand window
        self.ligShownEdit = self.findChild(QtWidgets.QLineEdit, 'titleLineEdit')
        self.bondingAtoms = self.findChild(QtWidgets.QLineEdit, 'bondingAtomLineEdit')
        self.ligWindow = self.findChild(gl.GLViewWidget,'LigandWidget')
        self.ligWindow = MyView(self.ligWindow)
        self.ligWindow.setGeometry(10,10,350,320)
        self.ligWindow.setCameraPosition(distance=10)
        self.addLigand = self.findChild(QtWidgets.QPushButton,'addNewLigandButton')
        self.addLigand.clicked.connect(self.addLigandButtonPushed)
        self.smilesLineEdit = self.findChild(QtWidgets.QLineEdit, 'smilesLineEdit')
        self.showSmilesBelowCheckBox = self.findChild(QtWidgets.QCheckBox, 'showSmilesBelowCheckBox')
        self.showSmilesBelowCheckBox.stateChanged.connect(self.showSmilesInLigandTab)

        #constraints tab
        self.chargeLineEdit = self.findChild(QtWidgets.QLineEdit, 'chargeLineEdit')
        self.uMultiplicity = self.findChild(QtWidgets.QLineEdit, 'uMultiplicity')
        #self.replaceMetalButton = self.findChild(QtWidgets.QCheckBox, 'checkBox')
        self.metalLabel = self.findChild(QtWidgets.QLabel, 'label_2')
        self.userMCEntry = self.findChild(QtWidgets.QLineEdit, 'user_mc')
        self.userCoordNum = self.findChild(QtWidgets.QLineEdit, 'userCoordNum')
        #frozen tab
        self.frozenTable = self.findChild(QtWidgets.QTableWidget,'tableWidget')
        self.addRowButton = self.findChild(QtWidgets.QPushButton, 'newFreezeRow')
        self.addRowButton.clicked.connect(self.addRow)
        self.removeRowButton = self.findChild(QtWidgets.QPushButton, 'removeFreezeRow')
        self.removeRowButton.clicked.connect(self._removeRow)
        
        #replace atom button, might replace this with something else
        self.replaceButton = self.findChild(QtWidgets.QPushButton, 'pushButton_2')
        self.replaceButton.clicked.connect(self.replaceAtoms)
        #next line doesn't work
        self.showAtomNumButton = self.findChild(QtWidgets.QLineEdit, 'actionViewatomic_num')
        
        #self.showAtomNumButton.triggered.connect(self.showAtomNum)
        global atom_type_edit
        global atomic_num
        global x_line_edit
        global z_line_edit
        global y_line_edit
        #self.selectMultipleButton = self.findChild(QtWidgets.QAction, 'actionSelect_Multiple')
        #self.selectMultipleButton.triggered.connect(self.selectMultiple)
        atom_type_edit = self.findChild(Qt.QLineEdit,'atomTypeLineEdit')
        atomic_num = self.findChild(Qt.QLineEdit,'atomNumberLineEdit')
        x_line_edit = self.findChild(Qt.QLineEdit,'xLineEdit')
        y_line_edit = self.findChild(Qt.QLineEdit,'yLineEdit')
        z_line_edit = self.findChild(Qt.QLineEdit,'zLineEdit')
        load_project_structure(self,'Libs', self.projectTreeWidget) #connects library with tree
        self.show()   

    def treeItemClicked(self, item, col):
        currFilePath = item.data(0, QtCore.Qt.UserRole) #gives us the path of the item clicked on
        currFile = os.path.basename(currFilePath) #removes the path
        currFile = os.path.splitext(currFile)[0] #removes the file extension
        self.ligShownEdit.setText(item.text(col)) 
        #delete everything in ligTabVisual so it doesn't get filled
        folder = 'ligTabVisual'
        for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete %s. Reason: %s' % (file_path, e))
        #direct stderr to a text file instead of null so you can see what went wrong if anything
        os.system('obabel ' + currFilePath +' -O ligTabVisual\\' + currFile + '.xyz -oxyz -m --gen3d  2>ligTabVisual\\obDialogue.txt')
        self.ligWindow.clear() #clears anything that might've been in the ligand window
        self.drawGraph('ligTabVisual\\'+currFile+'1.xyz', 'ligand') #hardcoding the 1 to make it work not sure why ob adds a 1


    def loadMolecule(self):
            self.openDialogBox()

    def openDialogBox(self):
        filename = QFileDialog.getOpenFileName()
        global template #so that the file can be accessed from any function
        template = filename[0]
        self.molWindow.clear()
        self.drawGraph(filename[0], 'main')

    def drawGraph(self, filename, mainOrLigand):
        import numpy as np
        iF = open(filename, 'r')
        iF.readline()
        iF.readline()
        line = iF.readline()
        if mainOrLigand=='main': #need to clear these everytime a new structure is drawn in the main window
            object_map.clear(); x.clear(); y.clear(); z.clear(); atom_list.clear(); atom_type.clear(); atom_Num.clear()
            while line:
                line = line.split()
                atom_type.append(line[0])
                atom_Num.append(ob.GetAtomicNum(line[0])) #gets atomnum from ob
                x.append(float(line[1]))
                y.append(float(line[2]))
                z.append(float(line[3]))
                line = iF.readline()
            iF.close()
            for i in range(len(x)):
                md = gl.MeshData.sphere(rows=10, cols = 20, radius = atom_radius_lib[str(atom_type[i])])
                #md = gl.MeshData.sphere(rows=10, cols = 20, radius = ob.GetCovalentRad(atom_Num[i]))
                m = gl.GLMeshItem(meshdata=md, smooth = False)
                m.setColor(QtGui.QColor(atom_type_to_color.get(atom_type[i],'pink'))) #color is pink if not in atom_type_to_color
                atom_color = atom_type_to_color.get(atom_type[i],'pink')
                #color = ob.vector3()
                #atom_color = ob.GetRGB(atom_Num[i], color.x(), color.y(), color.z())
                ##atom_color = ob.GetRGB(atom_Num[i], ob.vector3().x(), ob.vector3().y(), ob.vector3().z())
                #m.setColor(QtGui.QColor(atom_color))
                #m.setText(1)
                m.resetTransform()
                m.translate(x[i], y[i], z[i], local=True)
                #added an m element here
                #atom_list.append(Atom(atom_type[i],x[i],y[i],z[i], atom_radius_lib[str(atom_type[i])], m, atom_color))
                atom_list.append(Atom(atom_type[i],x[i],y[i],z[i], 0.2*ob.GetVdwRad(atom_Num[i]), m, atom_color))
                object_map[m] = i
                self.molWindow.addItem(m)
            calcLinesToDraw(self,x,y,z, atom_type, atom_Num, 'main')
        else:
            l_object_map.clear(); l_atom_list.clear(); l_atom_type.clear(); l_atom_Num.clear(); l_x.clear(); l_y.clear(); l_z.clear()
            while line:
                line = line.split()
                l_atom_type.append(line[0])
                l_atom_Num.append(ob.GetAtomicNum(line[0])) #gets atomnum from ob
                l_x.append(float(line[1]))
                l_y.append(float(line[2]))
                l_z.append(float(line[3]))
                line = iF.readline()
            iF.close()
            for i in range(len(l_x)):
                #md = gl.MeshData.sphere(rows=10, cols = 20, radius = atom_radius_lib[str(l_atom_type[i])])
                md = gl.MeshData.sphere(rows=10, cols = 20, radius = 0.2*ob.GetVdwRad(l_atom_Num[i]))
                m = gl.GLMeshItem(meshdata=md, smooth = False)
                m.setColor(QtGui.QColor(atom_type_to_color.get(l_atom_type[i],'pink'))) #color is pink if not in atom_type_to_color
                atom_color = atom_type_to_color.get(l_atom_type[i],'pink')
                #color = ob.vector3()
                #atom_color = ob.GetRGB(l_atom_Num[i], color.x(), color.y(), color.z())
                ##atom_color = ob.GetRGB(l_atom_Num[i], ob.vector3().x(), ob.vector3().y(), ob.vector3().z())
                #m.setColor(QtGui.QColor(atom_color))
                #m.setText(1)
                m.resetTransform()
                m.translate(l_x[i], l_y[i], l_z[i], local=True)
                #added an m element here
                #l_atom_list.append(Atom(l_atom_type[i],l_x[i],l_y[i],l_z[i], atom_radius_lib[str(l_atom_type[i])], m, atom_color))
                l_atom_list.append(Atom(l_atom_type[i],l_x[i],l_y[i],l_z[i], 0.2*ob.GetVdwRad(l_atom_Num[i]), m, atom_color))
                l_object_map[m] = i
                self.ligWindow.addItem(m)
            calcLinesToDraw(self,l_x,l_y,l_z, l_atom_type, l_atom_Num, 'ligand')


    def showAtomNum(self):
        pass
    
    def addRow(self):
        rowCount = self.tableWidget.rowCount()
        self.tableWidget.insertRow(rowCount )

    def _removeRow(self):
        if self.tableWidget.rowCount() >0:
            self.tableWidget.removeRow(self.tableWidget.rowCount()-1)

    def showSmilesInLigandTab(self, state):
        if state == QtCore.Qt.Checked:
            self.titleLineEdit.setText('') #clear the title of whatever was in it
            iF = open('ligTabVisual\\user.smi','w') #have to create a file for openbabel to read (I think??)
            iF.write(self.smilesLineEdit.text())
            iF.close()
            os.system('obabel ligTabVisual\\user.smi -O ligTabVisual\\user.xyz -oxyz -m --gen3d 2>ligTabVisual\\obDialogue.txt')
            self.ligWindow.clear() 
            self.drawGraph('ligTabVisual\\user1.xyz', 'ligand')

    def addLigandButtonPushed(self):
        msg = QMessageBox()
        message = ''
        if self.smilesLineEdit.text() == '' or self.ligShownEdit.text() == '' or self.bondingAtoms.text() == '':
            msg.setWindowTitle('please enter the missing info')
            if self.smilesLineEdit.text() == '': message = 'smiles format of the ligand\n'
            if self.ligShownEdit.text() == '' : message = message + 'title of the ligand\n'
            if self.bondingAtoms.text() == '': message = message + 'bonding atoms'
            msg.setText(message)
            x = msg.exec_()
        else:
            msg.setWindowTitle("is this correct?")
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            message = 'smiles format: '  + self.smilesLineEdit.text() + '\ntitle: ' + self.ligShownEdit.text() + '\nbonding atoms: ' + self.bondingAtoms.text()
            message = message + '\nafter clicking okay, your ligand can be found under the libraries tab in the user ligands folder'
            msg.setText(message)
            x = msg.exec_()
            if x == QMessageBox.Ok: 
                if not os.path.exists('Libs/user Ligands'): 
                    os.makedirs('Libs/user Ligands')
                userTitle = 'Libs/user Ligands/' + self.ligShownEdit.text()+ '.smi'
                f = open(userTitle,'w')
                f.write(self.smilesLineEdit.text()+ '\t' + self.bondingAtoms.text())
                self.projectTreeWidget.clear()
                load_project_structure(self,'Libs', self.projectTreeWidget) 


    def readFreezeTable(self):
        rowCount = self.tableWidget.rowCount()
        columnCount = self.tableWidget.columnCount()
        rowData = ''
        for row in range(rowCount):
            if(row >= 1): #skips the first row and adds a comma after every freeze row
                rowData = rowData + ','
            for column in range(columnCount):
                widgetItem = self.tableWidget.item(row,column)
                #if(widgetItem and widgetItem.text and column==1):
                if(widgetItem and widgetItem.text() and column>=1 and column<=4):
                    rowData = rowData + '-' + widgetItem.text()
                elif(widgetItem and widgetItem.text()):
                    rowData = rowData + widgetItem.text()
        return rowData

    def isEmpty(self, nCoordNum, nCoreNum, nCharge, nMultiplicity, nFrozen):
        if (not nCoordNum) or (not nCoreNum) or (not nCharge) or (not nMultiplicity) or (not nFrozen) or (not template):
            missingArgs = 'You are missing the following:\n'
            emptyMsg = QMessageBox()
            emptyMsg.setWindowTitle('Can\'t replace yet\n')
            emptyMsg.setIcon(QMessageBox.Critical)
            if (not nCoordNum): missingArgs = missingArgs + 'coord number\n'
            if (not nCoreNum): missingArgs = missingArgs + 'core number\n'
            if (not nCharge): missingArgs = missingArgs + 'charge\n'
            if (not nMultiplicity): missingArgs = missingArgs + 'multiplicity\n'
            if (not nFrozen): missingArgs = missingArgs + 'freezes\n'
            try: #because template is only set if you load a molecule
                template
            except NameError:
                missingArgs = missingArgs + 'template\n'
            emptyMsg.setText(missingArgs)
            x = emptyMsg.exec_()
            print(template)
            return True
        else: return False

    def ligandsSelected(self): #iterates through the ligands folder (whih is a tree or filewidget?)
        iterator = QtGui.QTreeWidgetItemIterator(self.projectTreeWidget, QtGui.QTreeWidgetItemIterator.Checked)
        ligandList = []
        ligPathList = []
        while iterator.value():
            item = iterator.value()
            itemPath = item.data(0, QtCore.Qt.UserRole)
            if os.path.isfile(itemPath): 
                ligPathList.append(itemPath) 
                ligandList.append(item.text(0))
            iterator += 1
        return ligandList, ligPathList

    def updateUserLigFolder(self,ligandList,ligPathList):
        #delete all ligands in userLigands (a folder separate from the ligands inside of smarties)
        import shutil
        folder = 'Smarties\\userLigands'
        for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete %s. Reason: %s' % (file_path, e))
        #get just the name of template without the path and file extension
        templateName = os.path.basename(template)
        templateName = os.path.splitext(templateName)[0]
        #if not os.path.exists('my_folder'): #rename it if it does exist??
        dirpath = 'Mason Files Generated By GUI\\' + templateName #make a folder for all the ligands the user selected
        if os.path.exists(dirpath):
            counter = 1
            while os.path.exists(dirpath):
                dirpath = 'Mason Files Generated By GUI\\' + templateName + " (" + str(counter) + ")"
                counter += 1
            os.makedirs(dirpath)
        else: os.makedirs(dirpath) 
        #for each ligand, copy it and paste it into the user folder
        for i in range(len(ligandList)):
            shutil.copyfile(ligPathList[i],dirpath + '\\' + ligandList[i]) #(source, destination)
        return templateName
    
    def createConfigFile(self, templateName,nCoordNum, nCoreNum, nCharge, nMultiplicity, nFrozen):
        newFileContent = ''
        iF = open('Smarties\config.ini','r')
        for line in iF:
            strippedLine = line.strip()
            newLine=strippedLine
            if 'coordNum' in strippedLine:
                newLine = strippedLine.replace(strippedLine, 'coordNum = ' + nCoordNum)
            elif 'Template' in strippedLine:
                newLine = strippedLine.replace(strippedLine, 'Template = ' + template)
            elif 'coreNum' in strippedLine:
                newLine = strippedLine.replace(strippedLine, 'coreNum = ' + nCoreNum)
            elif 'charge' in strippedLine:
                newLine = strippedLine.replace(strippedLine, 'charge = ' + nCharge)
            elif 'multiplicity' in strippedLine:
                newLine = strippedLine.replace(strippedLine, 'multipliciity = ' + nMultiplicity)
            elif 'freeze' in strippedLine:
                newLine = strippedLine.replace(strippedLine, 'freeze = ' + nFrozen)
            elif 'LigandLibDir' in strippedLine: #should probably do exact location like how template is exact loaction
                newLine = strippedLine.replace(strippedLine, 'LigandLibDir = Mason Files Generated By GUI/' + templateName + '/')
                #########might have to make this get the full path ^^^^^^
            newFileContent += newLine +'\n'
        iF.close()
        
        writingFile = open('Mason Files Generated By GUI\\'+ templateName + 'Config.ini', 'w')
        writingFile.write(newFileContent)
        writingFile.close()

    def replaceAtoms(self):
        nCoordNum = self.userCoordNum.text()
        nCoreNum = self.userMCEntry.text()
        nCharge = self.chargeLineEdit.text()
        nMultiplicity = self.uMultiplicity.text()
        nFrozen = self.readFreezeTable()
        
        if not self.isEmpty(nCoordNum, nCoreNum, nCharge, nMultiplicity, nFrozen): 
            msg = QMessageBox()
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            msg.setWindowTitle('Proceed to ligand replacement')
            msg.setIcon(QMessageBox.Question)
            msg.setText('Is this correct?\ncoordNum = ' + nCoordNum
            +'\nTemplate = ' + template + '\ncoreNum = ' + nCoreNum
            +'\ncharge = '+ nCharge + '\nmultiplicity = ' + nMultiplicity +'\nfreeze = ' + nFrozen)
            x = msg.exec_()

        if x == QMessageBox.Ok:
            ligandList, ligPathList = self.ligandsSelected() #get the ligands that the user selects
            templateName = self.updateUserLigFolder(ligandList, ligPathList)
            self.createConfigFile(templateName, nCoordNum, nCoreNum, nCharge, nMultiplicity, nFrozen)
              

def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    print(atom_type_edit)
    sys.exit(app.exec_())

if __name__ == '__main__':         
    main()
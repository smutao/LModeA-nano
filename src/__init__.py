'''
LModeA-nano: Local Vibrational Mode Analysis for Solids and Molecules
as a PyMOL Plugin

Ver 1.0.1 - Feb 26, 2022

Author: Yunwen Tao, Ph.D.

References:
1) LModeA-nano: A PyMOL Plugin for Calculating Bond Strength in Solids, Surfaces, and Molecules via Local Vibrational Mode Analysis, Y. Tao, W. Zou, S. Nanayakkara, and E. Kraka, J. Chem. Theory Comput., in press (2022)
2) In Situ Measure of Intrinsic Bond Strength in Crystalline Structures: Local Vibrational Mode Theory for Periodic Systems, Y. Tao, W. Zou, D. Sethio, N. Verma, Y. Qiu, C. Tian, D. Cremer, and E. Kraka, J. Chem. Theory Comput., 15, 1761-1776 (2019)

License: BSD-2-Clause
'''

from __future__ import absolute_import
from __future__ import print_function

# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.

# numpy is required for this plugin

import os
import os.path
from pymol.wizard import Wizard # import Wizard class 
from pymol import cmd


global_delocalized_bonds_list = []
global_natom = 0
global_elem = []
form = None
au2dy=15.56893
au2dy_a=4.359744
au2wn=5140.48715246
global_fm = [] # Hessian matrix
global_dim = 0
global_mass = []
global_minv = []
global_row_index=0
global_save_text = []
global_iflinmol = False



def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('LModeA-nano', run_plugin_gui)



from pymol.Qt import QtWidgets
from PyQt5.QtWidgets import QTableWidget,QTableWidgetItem # added by YTAO  
from PyQt5.QtCore import Qt # added     

#
object_prefix = "_pw"
class angleWizard(Wizard):

    global form
    def __init__(self):
        Wizard.__init__(self)
        
        self.pick_count = 0
        self.object_prefix = object_prefix
        
        self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
        cmd.set("mouse_selection_mode",0) # set selection mode to atomic
        cmd.deselect()

    def reset(self):
        cmd.delete(self.object_prefix + "*") # removed the selection of atoms "_pw*"
        #cmd.delete("sele*")
        cmd.delete("_indicate*")
        cmd.unpick()
        self.pick_count = 0

        cmd.refresh_wizard()

    def cleanup(self):
        cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode
        self.reset()
        #self.delete_all()

    def get_prompt(self): # print some info on the main screen
        self.prompt = None
        if self.pick_count == 0:
            self.prompt = [ 'Please click on the first angle atom...']
        elif self.pick_count == 1:
            self.prompt = [ 'Please click on the second angle atom...' ]
        elif self.pick_count == 2:
            self.prompt = [ 'Please click on the third angle atom...' ]
        elif self.pick_count == 3:
            self.prompt = [ 'Angle selected, please click "Done" button.' ]
        return self.prompt

    def do_select(self, name):
        # "edit" only this atom, and not others with the object prefix
            cmd.edit("%s and not %s*" % (name, self.object_prefix))
            self.do_pick(0)        

    def pickNextAtom(self, atom_name):
        # transfer the click selection to a named selection
        cmd.select(atom_name, "(pk1)")
        #print(atom_name) # _pw0, _pw1, _pw2 ...

        # delete the click selection
        cmd.unpick()

        # using the magic of indicate, highlight stuff
        indicate_selection = "_indicate" + self.object_prefix # starting with "_" can hide from window 

        #print(indicate_selection)
        cmd.select(indicate_selection, atom_name)
        cmd.enable(indicate_selection)

    def do_pick(self, picked_bond):

        # this shouldn't actually happen if going through the "do_select"
        if picked_bond:
            self.error = "Error: please select bonds, not atoms"
            print(self.error)
            return

        if self.pick_count>=3:
           self.pick_count=0

        atom_name = self.object_prefix + str(self.pick_count)

        if self.pick_count < 3:
            self.pickNextAtom(atom_name)
            #if self.pick_count == 0:
            #   self.pick_count += 1
            #   self.error = None
               # necessary to force update of the prompt
            #   cmd.refresh_wizard()
            self.pick_count += 1
            cmd.refresh_wizard()
        else:
            self.pickNextAtom(atom_name)  
            #self.pick_count = 0
            self.pick_count += 1
            cmd.refresh_wizard()
            self.reset()

    def done(self):
        global global_elem
        global global_natom
        global global_fm
        global global_dim
        global global_minv

        #print("an angle has been selected...")

        obj0 = cmd.get_object_list("_pw0")[0]
        obj1 = cmd.get_object_list("_pw1")[0]     
        obj2 = cmd.get_object_list("_pw2")[0]
      
        if obj0 == "geom" and obj1 == "supercell" and obj2 == "supercell":
           id0 = cmd.identify("_pw0")[0]
           id1 = cmd.identify("_pw1")[0]
           id1_map = id1%global_natom
           id2 = cmd.identify("_pw2")[0]
           id2_map = id2%global_natom

           car0=cmd.get_model('_pw0',1).get_coord_list()[0]
           car1=cmd.get_model('_pw1',1).get_coord_list()[0]
           car2=cmd.get_model('_pw2',1).get_coord_list()[0]
           a = round(calc_ang(car0,car1,car2),2)

           elem0 = global_elem[id0-1]
           elem1 = global_elem[id1_map-1]
           elem2 = global_elem[id2_map-1]
           if id1_map == 0:
              id1_map = global_natom 
           if id2_map == 0:
              id2_map = global_natom
           bvec = bmat_angle(car0,car1,car2,global_natom,id0,id1_map,id2_map)
           if id1_map == id0 or id2_map == id0:
              pass
           else:
              q_str = elem0+str(id0)+"-"+elem1+str(id1)+"("+str(id1_map)+")"+"-"+elem2+str(id2)+"("+str(id2_map)+")"
              rowPosition = form.tableWidget.rowCount()
              form.tableWidget.insertRow(rowPosition)
               
              #it=QTableWidgetItem(str(rowPosition+1))
              #it.setTextAlignment(Qt.AlignCenter)
              #it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              #form.tableWidget.setItem(rowPosition,0,it)

              it=QTableWidgetItem(q_str)
              it.setTextAlignment(Qt.AlignCenter)
              it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,0,it)

              it_1 = QTableWidgetItem(str(a))
              it_1.setTextAlignment(Qt.AlignCenter)
              it_1.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,1,it_1)

              ka,freq = calc_fc(global_fm,global_minv,bvec,int(global_dim))
              ka = round(ka * au2dy_a,3)
              it_2 = QTableWidgetItem(str(ka))
              it_2.setTextAlignment(Qt.AlignCenter)
              it_2.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,2,it_2)

              freq = int(round(freq*au2wn))
              #print("freq:",freq)
              it_3 = QTableWidgetItem(str(freq))
              it_3.setTextAlignment(Qt.AlignCenter)
              it_3.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,3,it_3)



        if obj0 == "geom" and obj1 == "supercell" and obj2 == "geom":
           id0 = cmd.identify("_pw0")[0]
           id1 = cmd.identify("_pw1")[0]
           id1_map = id1%global_natom
           id2 = cmd.identify("_pw2")[0]
 
           car0=cmd.get_model('_pw0',1).get_coord_list()[0]
           car1=cmd.get_model('_pw1',1).get_coord_list()[0]
           car2=cmd.get_model('_pw2',1).get_coord_list()[0]
           a = round(calc_ang(car0,car1,car2),2)

           elem0 = global_elem[id0-1]
           elem1 = global_elem[id1_map-1]
           elem2 = global_elem[id2-1]
           if id1_map == 0:
              id1_map = global_natom 
           bvec = bmat_angle(car0,car1,car2,global_natom,id0,id1_map,id2)
           if id1_map == id0 or id1_map == id2:
              pass
           else:
              q_str = elem0+str(id0)+"-"+elem1+str(id1)+"("+str(id1_map)+")"+"-"+elem2+str(id2)
              rowPosition = form.tableWidget.rowCount()
              form.tableWidget.insertRow(rowPosition)
               
              #it=QTableWidgetItem(str(rowPosition+1))
              #it.setTextAlignment(Qt.AlignCenter)
              #it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              #form.tableWidget.setItem(rowPosition,0,it)

              it=QTableWidgetItem(q_str)
              it.setTextAlignment(Qt.AlignCenter)
              it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,0,it)

              it_1 = QTableWidgetItem(str(a))
              it_1.setTextAlignment(Qt.AlignCenter)
              it_1.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,1,it_1)

              ka,freq = calc_fc(global_fm,global_minv,bvec,int(global_dim))
              ka = round(ka * au2dy_a,3)
              it_2 = QTableWidgetItem(str(ka))
              it_2.setTextAlignment(Qt.AlignCenter)
              it_2.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,2,it_2)


              freq = int(round(freq*au2wn))
              it_3 = QTableWidgetItem(str(freq))
              it_3.setTextAlignment(Qt.AlignCenter)
              it_3.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,3,it_3)



        if obj0 == "geom" and obj1 == "geom" and obj2 == "supercell":
           id0 = cmd.identify("_pw0")[0]
           id1 = cmd.identify("_pw1")[0]
           id2 = cmd.identify("_pw2")[0]
           id2_map = id2%global_natom

           car0=cmd.get_model('_pw0',1).get_coord_list()[0]
           car1=cmd.get_model('_pw1',1).get_coord_list()[0]
           car2=cmd.get_model('_pw2',1).get_coord_list()[0]
           a = round(calc_ang(car0,car1,car2),2)
           
           elem0 = global_elem[id0-1]
           elem1 = global_elem[id1-1]
           elem2 = global_elem[id2_map-1]
           if id2_map == 0:
              id2_map = global_natom
           bvec = bmat_angle(car0,car1,car2,global_natom,id0,id1,id2_map) 

           if id2_map == id0 or id2_map == id1:
              pass
           else:
              q_str = elem0+str(id0)+"-"+elem1+str(id1)+"-"+elem2+str(id2)+"("+str(id2_map)+")"  
              rowPosition = form.tableWidget.rowCount()
              form.tableWidget.insertRow(rowPosition)
               
              #it=QTableWidgetItem(str(rowPosition+1))
              #it.setTextAlignment(Qt.AlignCenter)
              #it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              #form.tableWidget.setItem(rowPosition,0,it)

              it=QTableWidgetItem(q_str)
              it.setTextAlignment(Qt.AlignCenter)
              it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,0,it)

              it_1 = QTableWidgetItem(str(a))
              it_1.setTextAlignment(Qt.AlignCenter)
              it_1.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,1,it_1)

              ka,freq = calc_fc(global_fm,global_minv,bvec,int(global_dim))
              ka = round(ka * au2dy_a,3)
              it_2 = QTableWidgetItem(str(ka))
              it_2.setTextAlignment(Qt.AlignCenter)
              it_2.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,2,it_2)


              freq = int(round(freq*au2wn))
              #print("freq:",freq)
              it_3 = QTableWidgetItem(str(freq))
              it_3.setTextAlignment(Qt.AlignCenter)
              it_3.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,3,it_3)


        if obj0 == "geom" and obj1 == "geom" and obj2 == "geom":
           id0 = cmd.identify("_pw0")[0]
           id1 = cmd.identify("_pw1")[0]
           id2 = cmd.identify("_pw2")[0]

           car0=cmd.get_model('_pw0',1).get_coord_list()[0]
           car1=cmd.get_model('_pw1',1).get_coord_list()[0]
           car2=cmd.get_model('_pw2',1).get_coord_list()[0]
           #r = round(calc_dis(car0,car1),4)
           a = round(calc_ang(car0,car1,car2),2)
           #print(r)
           elem0 = global_elem[id0-1]
           elem1 = global_elem[id1-1]
           elem2 = global_elem[id2-1]

           #bvec = bmat_bond(car0,car1,global_natom,id0,id1)
           #print(bvec)
           bvec = bmat_angle(car0,car1,car2,global_natom,id0,id1,id2)
           #print(bvec)

   
         
           q_str = elem0+str(id0)+"-"+elem1+str(id1)+"-"+elem2+str(id2)   #C1-O2 
           #print(q_str) 
           rowPosition = form.tableWidget.rowCount()
           #print(rowPosition) 
           form.tableWidget.insertRow(rowPosition)
           
           #it=QTableWidgetItem(str(rowPosition+1))
           #it.setTextAlignment(Qt.AlignCenter)
           #it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
           #form.tableWidget.setItem(rowPosition,0,it)

           it=QTableWidgetItem(q_str)
           it.setTextAlignment(Qt.AlignCenter)
           it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
           form.tableWidget.setItem(rowPosition,0,it)

           it_1 = QTableWidgetItem(str(a))
           it_1.setTextAlignment(Qt.AlignCenter)
           it_1.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
           form.tableWidget.setItem(rowPosition,1,it_1)
 
           ka,freq = calc_fc(global_fm,global_minv,bvec,int(global_dim))
           ka = round(ka * au2dy_a,3)
           it_2 = QTableWidgetItem(str(ka))
           it_2.setTextAlignment(Qt.AlignCenter)
           it_2.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
           form.tableWidget.setItem(rowPosition,2,it_2)


           freq = int(round(freq*au2wn))
           #print("freq:",freq)
           it_3 = QTableWidgetItem(str(freq))
           it_3.setTextAlignment(Qt.AlignCenter)
           it_3.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
           form.tableWidget.setItem(rowPosition,3,it_3)


        cmd.set_wizard()

    def get_panel(self): # show menu 
        return [
            [ 1, 'Angle Selection Wizard',''],
            [ 2, 'Reset','cmd.get_wizard().reset()'],
            [ 2, 'Done','cmd.get_wizard().done()'],
            [ 2, 'Close','cmd.set_wizard()'],
        ]



class bondWizard(Wizard): # class inheritance
    #from pymol.Qt import QtWidgets
    #from PyQt5.QtWidgets import QTableWidget,QTableWidgetItem # added by YTAO  
    #from PyQt5.QtCore import Qt # added     



    global form


    def __init__(self):
        Wizard.__init__(self)
        
        self.pick_count = 0
        self.object_prefix = object_prefix
        
        self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
        cmd.set("mouse_selection_mode",0) # set selection mode to atomic
        cmd.deselect()


    def reset(self):
        cmd.delete(self.object_prefix + "*") # removed the selection of atoms "_pw*"
        #cmd.delete("sele*")
        cmd.delete("_indicate*")
        cmd.unpick()
        self.pick_count = 0

        cmd.refresh_wizard()

    def cleanup(self):
        cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode
        self.reset()
        #self.delete_all()

    def get_prompt(self): # print some info on the main screen
        self.prompt = None
        if self.pick_count == 0:
            self.prompt = [ 'Please click on the first bond atom...']
        elif self.pick_count == 1:
            self.prompt = [ 'Please click on the second bond atom...' ]
        elif self.pick_count == 2:
            self.prompt = [ 'Bond selected, please click "Done" button.' ]
        return self.prompt


    def do_select(self, name):
        # "edit" only this atom, and not others with the object prefix
            cmd.edit("%s and not %s*" % (name, self.object_prefix))
            self.do_pick(0)        

    def pickNextAtom(self, atom_name):
        # transfer the click selection to a named selection
        cmd.select(atom_name, "(pk1)")
        #print(atom_name) # _pw0, _pw1, ...

        # delete the click selection
        cmd.unpick()

        # using the magic of indicate, highlight stuff
        indicate_selection = "_indicate" + self.object_prefix # starting with "_" can hide from window 

        #print(indicate_selection)
        cmd.select(indicate_selection, atom_name)
        cmd.enable(indicate_selection)

    def do_pick(self, picked_bond):

        # this shouldn't actually happen if going through the "do_select"
        if picked_bond:
            self.error = "Error: please select bonds, not atoms"
            print(self.error)
            return

        if self.pick_count>=2:
           self.pick_count=0

        atom_name = self.object_prefix + str(self.pick_count)

        if self.pick_count < 2:
            self.pickNextAtom(atom_name)
            #if self.pick_count == 0:
            #   self.pick_count += 1
            #   self.error = None
               # necessary to force update of the prompt
            #   cmd.refresh_wizard()
            self.pick_count += 1
            cmd.refresh_wizard()
        else:
            self.pickNextAtom(atom_name)  
            #self.pick_count = 0
            self.pick_count += 1
            cmd.refresh_wizard()
            self.reset()

    def done(self):
        global global_elem
        global global_natom
        global global_fm
        global global_dim
        global global_minv

        #print("a bond has been selected...")

        obj0 = cmd.get_object_list("_pw0")[0]
        obj1 = cmd.get_object_list("_pw1")[0]     
        if obj1 == "geom" and obj0 == "supercell":
           id0 = cmd.identify("_pw1")[0] # somehow
           id1 = cmd.identify("_pw0")[0] # confusing
           id1_map = id1%global_natom

           car1=cmd.get_model('_pw0',1).get_coord_list()[0] # id0-car0, id1-car1, 
           car0=cmd.get_model('_pw1',1).get_coord_list()[0] # easy to calculate B 
           r = round(calc_dis(car0,car1),4)

           elem0 = global_elem[id0-1]
           elem1 = global_elem[id1_map-1]
           if id1_map == 0:
              id1_map = global_natom

           bvec = bmat_bond(car0,car1,global_natom,id0,id1_map)
           #print(bvec)

           if id0 == id1_map:
              pass
           else:
              q_str = elem0+str(id0)+"-"+elem1+str(id1)+"("+str(id1_map)+")"
              rowPosition = form.tableWidget.rowCount()
              form.tableWidget.insertRow(rowPosition)
           
              #it=QTableWidgetItem(str(rowPosition+1))
              #it.setTextAlignment(Qt.AlignCenter)
              #it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              #form.tableWidget.setItem(rowPosition,0,it)

              it=QTableWidgetItem(q_str)
              it.setTextAlignment(Qt.AlignCenter)
              it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,0,it)

              it_1 = QTableWidgetItem(str(r))
              it_1.setTextAlignment(Qt.AlignCenter)
              it_1.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,1,it_1)
 
              ka,freq = calc_fc(global_fm,global_minv,bvec,int(global_dim)) 
              ka = round(ka * au2dy,3) # convert from atomic unit 
              it_2 = QTableWidgetItem(str(ka))
              it_2.setTextAlignment(Qt.AlignCenter)
              it_2.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,2,it_2)

              freq = int(round(freq*au2wn))
              #print("freq:",freq)
              it_3 = QTableWidgetItem(str(freq))
              it_3.setTextAlignment(Qt.AlignCenter)
              it_3.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,3,it_3)


        if obj0 == "supercell" and obj1 == "supercell":
           id0 = cmd.identify("_pw0")[0]
           id0_map = id0%global_natom
           id1 = cmd.identify("_pw1")[0]
           id1_map = id1%global_natom
           
           car0=cmd.get_model('_pw0',1).get_coord_list()[0]
           car1=cmd.get_model('_pw1',1).get_coord_list()[0]
           r = round(calc_dis(car0,car1),4)
            
           elem0 = global_elem[id0_map-1]
           elem1 = global_elem[id1_map-1]
           if id0_map == 0:
              id0_map = global_natom 
           if id1_map == 0:
              id1_map = global_natom 

           bvec = bmat_bond(car0,car1,global_natom,id0_map,id1_map)
           #print(bvec)

           q_str = elem0+str(id0_map)+"-"+elem1+str(id1_map)
           rowPosition = form.tableWidget.rowCount()
           form.tableWidget.insertRow(rowPosition)
           
           #it=QTableWidgetItem(str(rowPosition+1))
           #it.setTextAlignment(Qt.AlignCenter)
           #it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
           #form.tableWidget.setItem(rowPosition,0,it)

           it=QTableWidgetItem(q_str)
           it.setTextAlignment(Qt.AlignCenter)
           it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
           form.tableWidget.setItem(rowPosition,0,it)

           it_1 = QTableWidgetItem(str(r))
           it_1.setTextAlignment(Qt.AlignCenter)
           it_1.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
           form.tableWidget.setItem(rowPosition,1,it_1)
 
           ka,freq = calc_fc(global_fm,global_minv,bvec,int(global_dim)) 
           ka = round(ka * au2dy,3) # convert from atomic unit 
           it_2 = QTableWidgetItem(str(ka))
           it_2.setTextAlignment(Qt.AlignCenter)
           it_2.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
           form.tableWidget.setItem(rowPosition,2,it_2)


           freq = int(round(freq*au2wn))
           #print("freq:",freq)
           it_3 = QTableWidgetItem(str(freq))
           it_3.setTextAlignment(Qt.AlignCenter)
           it_3.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
           form.tableWidget.setItem(rowPosition,3,it_3)




        if obj0 == "geom" and obj1 == "supercell":
           id0 = cmd.identify("_pw0")[0]
           id1 = cmd.identify("_pw1")[0]
           id1_map = id1%global_natom

           car0=cmd.get_model('_pw0',1).get_coord_list()[0]
           car1=cmd.get_model('_pw1',1).get_coord_list()[0]
           r = round(calc_dis(car0,car1),4)
           
           elem0 = global_elem[id0-1]
           elem1 = global_elem[id1_map-1]
           if id1_map == 0:
              id1_map = global_natom 
           bvec = bmat_bond(car0,car1,global_natom,id0,id1_map)
           #print(bvec)

           if id1_map == id0:
              pass
           else:
              q_str = elem0+str(id0)+"-"+elem1+str(id1)+"("+str(id1_map)+")"# C1-O12(6)
              rowPosition = form.tableWidget.rowCount()
              form.tableWidget.insertRow(rowPosition)
           
              #it=QTableWidgetItem(str(rowPosition+1))
              #it.setTextAlignment(Qt.AlignCenter)
              #it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              #form.tableWidget.setItem(rowPosition,0,it)

              it=QTableWidgetItem(q_str)
              it.setTextAlignment(Qt.AlignCenter)
              it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,0,it)

              it_1 = QTableWidgetItem(str(r))
              it_1.setTextAlignment(Qt.AlignCenter)
              it_1.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,1,it_1)

              ka,freq = calc_fc(global_fm,global_minv,bvec,int(global_dim)) 
              ka = round(ka * au2dy,3) # convert from atomic unit 
              it_2 = QTableWidgetItem(str(ka))
              it_2.setTextAlignment(Qt.AlignCenter)
              it_2.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,2,it_2)


              freq = int(round(freq*au2wn))
              #print("freq:",freq)
              it_3 = QTableWidgetItem(str(freq))
              it_3.setTextAlignment(Qt.AlignCenter)
              it_3.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
              form.tableWidget.setItem(rowPosition,3,it_3)


        #print(obj0,obj1)
        if obj0 == "geom" and obj1 == "geom":
           id0 = cmd.identify("_pw0")[0]
           id1 = cmd.identify("_pw1")[0]
           car0=cmd.get_model('_pw0',1).get_coord_list()[0]
           car1=cmd.get_model('_pw1',1).get_coord_list()[0]
           r = round(calc_dis(car0,car1),4)
           #print(r)
           elem0 = global_elem[id0-1]
           elem1 = global_elem[id1-1]
           bvec = bmat_bond(car0,car1,global_natom,id0,id1)
           #print(bvec)
            
        
           q_str = elem0+str(id0)+"-"+elem1+str(id1)   #C1-O2 
           #print(q_str) 
           rowPosition = form.tableWidget.rowCount()
           #print(rowPosition) 
           form.tableWidget.insertRow(rowPosition)
           
           #it=QTableWidgetItem(str(rowPosition+1))
           #it.setTextAlignment(Qt.AlignCenter)
           #it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
           #form.tableWidget.setItem(rowPosition,0,it)

           it=QTableWidgetItem(q_str)
           it.setTextAlignment(Qt.AlignCenter)
           it.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
           form.tableWidget.setItem(rowPosition,0,it)

           it_1 = QTableWidgetItem(str(r))
           it_1.setTextAlignment(Qt.AlignCenter)
           it_1.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
           form.tableWidget.setItem(rowPosition,1,it_1)
           
           # then calculate Wilson B matrix and force constant 
           ka,freq = calc_fc(global_fm,global_minv,bvec,int(global_dim)) 
           ka = round(ka * au2dy,3) # convert from atomic unit 
           it_2 = QTableWidgetItem(str(ka))
           it_2.setTextAlignment(Qt.AlignCenter)
           it_2.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
           form.tableWidget.setItem(rowPosition,2,it_2)


           freq = int(round(freq*au2wn))
           #print("freq:",freq)
           it_3 = QTableWidgetItem(str(freq))
           it_3.setTextAlignment(Qt.AlignCenter)
           it_3.setFlags(  Qt.ItemIsSelectable |  Qt.ItemIsEnabled  )
           form.tableWidget.setItem(rowPosition,3,it_3)




        cmd.set_wizard()

    def get_panel(self): # show menu 
        return [
            [ 1, 'Bond Selection Wizard',''],
            [ 2, 'Reset','cmd.get_wizard().reset()'],
            #[ 2, 'Arm Atoms Selection Done', 'cmd.get_wizard().finish_1arm()'],
            #[ 2, 'Done','cmd.set_wizard()'],
            [ 2, 'Done','cmd.get_wizard().done()'],
            [ 2, 'Close','cmd.set_wizard()'],
        ]

# global functions 

def calc_minv(mass):
    natom = len(mass)
    minv = []

    for i in range(3*natom):
        minv.append([])
    for i in range(3*natom):
        for j in range(3*natom):
            minv[i].append(0)
    for i in range(natom):
        for j in range(3):
            m = 3*i+j
            minv[m][m] = 1.0/mass[i]

    return minv


def get_symbol(z):
     
     s = ["H","He",\
             "Li","Be","B","C","N","O","F","Ne",\
             "Na","Mg","Al","Si","P","S","Cl","Ar",\
             "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",\
             "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",\
             "Cs","Ba",     "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",\
                      "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",\
             "Fr","Ra",     "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",\
             "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"\
             ]

     return s[z-1]



def get_mass(elem):
    #average atomic mass
     mass=[ 1.007940000,   4.002602000,   6.941000000,
     9.012183100,  10.811000000,  12.010700000,
     14.006700000,  15.999400000,  18.998403163,
     20.179700000,  22.989769280,  24.305000000,
     26.981538500,  28.085500000,  30.973761998,
     32.065000000,  35.453000000,  39.948000000,
     39.098300000,  40.078000000,  44.955908000,
     47.867000000,  50.941500000,  51.996100000,
     54.938044000,  55.845000000,  58.933194000,
     58.693400000,  63.546000000,  65.380000000,
     69.723000000,  72.640000000,  74.921595000,
     78.971000000,  79.904000000,  83.798000000,
     85.467800000,  87.620000000,  88.905840000,
     91.224000000,  92.906370000,  95.950000000,
     98.907200000, 101.070000000, 102.905500000,
     106.420000000, 107.868200000, 112.414000000,
     114.818000000, 118.710000000, 121.760000000,
     127.600000000, 126.904470000, 131.293000000,
     132.905451960, 137.327000000, 138.905470000,
     140.116000000, 140.907660000, 144.242000000,
     144.900000000, 150.360000000, 151.964000000,
     157.250000000, 158.925350000, 162.500000000,
     164.930330000, 167.259000000, 168.934220000,
     173.054000000, 174.966800000, 178.490000000,
     180.947880000, 183.840000000, 186.207000000,
     190.230000000, 192.217000000, 195.084000000,
     196.966569000, 200.590000000, 204.383300000,
     207.200000000, 208.980400000, 208.982400000,
     209.987100000, 222.017600000, 223.019700000,
     226.024500000, 227.027700000, 232.037700000,
     231.035880000, 238.028910000, 237.048200000,
     239.064200000, 243.061400000, 247.070400000,
     247.070300000, 251.079600000, 252.083000000,
     257.059100000, 258.098400000, 259.101000000,
     262.109700000, 261.108800000, 262.114100000,
     266.121900000, 264.120100000, 265.000000000,
     268.138800000, 269.000000000, 272.000000000,
     277.000000000, 280.000000000, 280.000000000,
     280.000000000, 280.000000000, 280.000000000,
     280.000000000, 280.000000000, 280.000000000]

     s = ["H","He",\
             "Li","Be","B","C","N","O","F","Ne",\
             "Na","Mg","Al","Si","P","S","Cl","Ar",\
             "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",\
             "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",\
             "Cs","Ba",     "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",\
                      "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",\
             "Fr","Ra",     "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",\
             "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"\
             ]
     s_upper = [str1.upper() for str1 in s]
     #index = s.index(elem)
     index = s_upper.index(elem.upper())

     return mass[index]




def calc_fc(fm,minv,bv,dim):
    import numpy as np
    import numpy.linalg as linalg
    from numpy.linalg import inv,pinv
   
    global global_iflinmol

    global global_save_text
    to_save = []
    ifsave = 1
    if len(global_save_text) >= 7 :
       ifsave = 0

    ka = 0.0
    natom3 = len(fm)
    f = np.array(fm)

    minv = np.array(minv)

     
    if dim !=0: # solids
      val,vec = linalg.eig(f)
      idx = val.argsort()
      val=val[idx]
      val = val.real
      vec=vec[:,idx]
      

    else: # molecules
      #something goes wrong with the formula (17) from JCTC_2019 paper, therefore switch to (21)
      b=np.array(bv) 
      #print("b vec",b) 
      #ka = 1.0/(np.matmul(np.matmul(b,pinv(f)),np.transpose(b))) 
      val,vec = linalg.eig(f)
      idx = val.argsort()
      val=val[idx]
      val = val.real
      vec=vec[:,idx]
 
    for i in range(6):
       to_save.append( round(float(val[i]),8) )
    if ifsave == 1:
       global_save_text.append(to_save)




  

    #print("dim is",dim)
    if dim == 2 or dim == 3:
       L = vec[:,3:natom3]
    if dim == 1:
       L = vec[:,4:natom3]
    if dim == 0:
       #print("if lin mol:",global_iflinmol) 
       if global_iflinmol: 
          L = vec[:,5:natom3]
       else:
          L = vec[:,6:natom3]  
    
    if dim != 0 or dim == 0:
       K=np.matmul(np.transpose(L),np.matmul(fm,L))
       b=np.array(bv)
       D=np.matmul(b,L)
       #print("D is",D)
       ka = 1.0/(np.matmul(np.matmul(D,inv(K)),np.transpose(D)))
    
    # frequency
    g = np.matmul(np.matmul(b,minv),np.transpose(b))
    adia_m = 1.0/g
    tmp = ka/adia_m
    if tmp < 0:
       freq = -abs(tmp)**0.5 
    else:
       freq =  tmp**0.5 

    return ka.real,freq


def calc_ang(l1,l2,l3):
    from numpy import sign
    from math import acos

    t1 = list_minus(l1,l2)
    t2 = list_minus(l3,l2)
    d12 = list_dotp(t1,t1)**0.5
    d32 = list_dotp(t2,t2)**0.5
    t1 = list_scale(t1,1.0/d12)
    t2 = list_scale(t2,1.0/d32)
    cphi = list_dotp(t1,t2)
    if abs(cphi)>1.0:
       cphi = sign(cphi)*1.0
    angle=acos(cphi)*57.295779513


    return angle




def calc_dis(l1,l2):
    d = 0.0
    for i in range(3):
        d = d + (l1[i]-l2[i])**2
    d = d**0.5
    return d 


def list_minus(a,b): # a-b
    c = []
    if len(a) != len(b):
        return c
    for i in range(len(a)):
        c.append(a[i]-b[i])
    return c

def list_dotp(a,b):
    n = len(a)
    s =  0.0
    for i in range(n):
        s = s + a[i]*b[i]
    return s 
 
def list_add(a,b):
    c = []
    if len(a) != len(b):
        return c
    for i in range(len(a)):
        c.append(a[i]+b[i])
    return c

def list_scale(l1,s):
    l2 = []
    for i in range(len(l1)):
        l2.append(l1[i]*s)
    return l2

def bmat_angle(c1,c2,c3,natom,id1,id2,id3):
    ang2bohr = 1.8897259886
    b = []
    for i in range(natom):
        for j in range(3):
            b.append(0.0)
    
    t1 = list_minus(c1,c2)
    t2 = list_minus(c3,c2)
    d12 = list_dotp(t1,t1)**0.5
    d32 = list_dotp(t2,t2)**0.5
    t1 = list_scale(t1,1.0/d12)
    t2 = list_scale(t2,1.0/d32)
    cphi = list_dotp(t1,t2)
    sphi = (1.0-cphi*cphi)**0.5

    b1 = list_scale(t2,-1.0) 
    b1 = AccAB(cphi,t1,b1)
    cf = 1.0/(d12*sphi)
    b1 = list_scale(b1,cf)

    b3 = list_scale(t1,-1.0)
    b3 = AccAB(cphi,t2,b3)
    cf = 1.0/(d32*sphi)
    b3 = list_scale(b3,cf)

    b2 = list_scale(b1,-1.0)
    b2 = list_minus(b2,b3) 
   
    b[(id1-1)*3  ] = b1[0]/ang2bohr
    b[(id1-1)*3+1] = b1[1]/ang2bohr
    b[(id1-1)*3+2] = b1[2]/ang2bohr

    b[(id2-1)*3  ] = b2[0]/ang2bohr
    b[(id2-1)*3+1] = b2[1]/ang2bohr
    b[(id2-1)*3+2] = b2[2]/ang2bohr

    b[(id3-1)*3  ] = b3[0]/ang2bohr
    b[(id3-1)*3+1] = b3[1]/ang2bohr
    b[(id3-1)*3+2] = b3[2]/ang2bohr

    return b 


def AccAB(c,A,B):
    n=len(A)
    for i in range(n):
        B[i] = B[i] + c*A[i]

    return B   




def bmat_bond(c1,c2,natom,id1,id2):
    # id1 id2 start from 1
    b = []
    for i in range(natom):
        for j in range(3):
            b.append(0.0)
    r = calc_dis(c1,c2)
    diff = list_minus(c2,c1)

    b1 = list_scale(diff,-1.0/r)
    b2 = list_scale(diff, 1.0/r)

    b[(id1-1)*3  ] = b1[0]
    b[(id1-1)*3+1] = b1[1]
    b[(id1-1)*3+2] = b1[2]

    b[(id2-1)*3  ] = b2[0]
    b[(id2-1)*3+1] = b2[1]
    b[(id2-1)*3+2] = b2[2]

    return b



# global reference to avoid garbage collection of our dialog
dialog = None


def run_plugin_gui():
    '''
    Open our custom dialog
    '''
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()


def make_dialog():
    global form 


    # entry point to PyMOL's API
    from pymol import cmd

    # pymol.Qt provides the PyQt5 interface, but may support PyQt4
    # and/or PySide as well
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi
    from pymol.Qt.utils import getSaveFileNameWithExt
    from PyQt5.QtWidgets import QTableWidget,QTableWidgetItem # added by YTAO  
    from PyQt5.QtCore import Qt # added     


    #bondwiz = bondWizard()





    # create a new Window
    dialog = QtWidgets.QDialog()

    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'gui-2.ui')
    form = loadUi(uifile, dialog)

    bondwiz = bondWizard()
    anglewiz = angleWizard()



    # callback for the "Browse" button
    def browse_filename():
        filename = getSaveFileNameWithExt(
            dialog, 'Save As...', filter='Text File (*.txt)')
        if filename:
            #form.input_filename.setText(filename)
           return filename


    def getOpenFileNameWithExt(*args, **kwargs):
      """
      Return a file name, append extension from filter if no extension provided.
      """
      import os, re

      fname, filter = QtWidgets.QFileDialog.getOpenFileName(*args, **kwargs)

      if not fname:
          return ''

      if '.' not in os.path.split(fname)[-1]:
          m = re.search(r'\*(\.[\w\.]+)', filter)
          if m:
              # append first extension from filter
              fname += m.group(1)

      return fname


    def find_input_file():
        input_file_name = getOpenFileNameWithExt( 
              dialog, 'Open...',filter='Any Data File (*.*)')
        if input_file_name:
           form.input_geomname.setText(input_file_name)

        return

    def read_vasp(f1):
        global global_save_text


        posfile = ""
        fcfile = ""

        to_save=[]
 
        k=0
        with open(f1) as f:
            for line in f:
                if "@VASP" in line:
                    k = 1
                    to_save.append(line.strip())
        f.close()
        if k == 0:
           return posfile,fcfile
        with open(f1) as f:
            for line in f:
                if "POSCON" in line and "=" in line:
                    posfile = line.split("=")[-1].strip()
                    to_save.append(line.strip())
                if "FCPHONOPY" in line and "=" in line:
                    fcfile = line.split("=")[-1].strip()
                    to_save.append(line.strip())

        global_save_text.append(to_save)

        
        return posfile,fcfile


    def read_qe(f1):
        global global_save_text

        to_save = []

        mode = ""
        inpfile = ""
        fcfile = ""
        k=0
        with open(f1) as f:
            for line in f:
                if "@QE" in line:
                    k = 1
                    to_save.append(line.strip())
        f.close()
        if k == 0 :
           return mode,inpfile,fcfile
        with open(f1) as f:
            for line in f:
                if "mode" in line and "=" in line:
                    mode = line.split("=")[-1].strip()
                    to_save.append(line.strip())
                if "GeomIn" in line and "=" in line:
                    inpfile = line.split("=")[-1].strip()
                    to_save.append(line.strip())
                if "q2rFC" in line and "=" in line:
                    fcfile = line.split("=")[-1].strip()
                    to_save.append(line.strip())
                if "FCPHONOPY" in line and "=" in line:
                    fcfile = line.split("=")[-1].strip()
                    to_save.append(line.strip())
        
        global_save_text.append(to_save)

        return  mode,inpfile,fcfile

    
    def read_unimovib(f1):
        global global_save_text
        to_save = []

        umvfile = ""
        k=0
        with open(f1) as f:
            for line in f:
                if "@unimovib" in line:
                    k = 1
                    to_save.append(line.strip())

        f.close()
        if k == 0 :
           return umvfile
        with open(f1) as f:
            for line in f:
                if "UMV" in line and "=" in line:
                    umvfile = line.split("=")[-1].strip()
                    to_save.append(line.strip())
 
        global_save_text.append(to_save)
        return umvfile
    
    def read_castep(f1):
        global global_save_text
        to_save = []
        
        outfile = ""
        k=0
        with open(f1) as f:
            for line in f:
                if "@castep" in line:
                    k = 1
                    to_save.append(line.strip())

        f.close()
        #print("k=",k)
        if k == 0 :
           return outfile
        with open(f1) as f:
            for line in f:
                if "output" in line and "=" in line:
                    outfile = line.split("=")[-1].strip()
                    to_save.append(line.strip())
        global_save_text.append(to_save)
  
        return outfile 


    def read_cp2k(f1):
        global global_save_text
        to_save = []

        fcfile = ""
        inpfile = ""
        mode = ""
        
        k=0
        with open(f1) as f:
            for line in f:
                if "@cp2k" in line:
                    k = 1
                    to_save.append(line.strip())

        f.close()
        if k == 0 :
           return fcfile,inpfile
        with open(f1) as f:
            for line in f:
                if "mode" in line and "=" in line: # mode = native,phonopy
                    mode = line.split("=")[-1].strip()
                    to_save.append(line.strip())
                if "OUTPUT" in line and "=" in line:
                    fcfile = line.split("=")[-1].strip()
                    to_save.append(line.strip())
                if "INPUT" in line and "=" in line:
                    inpfile = line.split("=")[-1].strip()
                    to_save.append(line.strip())
                if "FCPHONOPY" in line and "=" in line:
                    fcfile = line.split("=")[-1].strip()
                    to_save.append(line.strip())

        global_save_text.append(to_save)

        return mode,fcfile,inpfile

    
    def read_crystal(f1):
        global global_save_text
        outfile = ""
        fcfile = ""
        
        to_save = []

        k=0
        with open(f1) as f:
            for line in f:
                if "@crystal" in line:
                    k = 1
                    to_save.append(line.strip())
        f.close()
        if k == 0:
           return outfile,fcfile
        with open(f1) as f:
            for line in f:
                if "OUTPUT" in line and "=" in line:
                    outfile = line.split("=")[-1].strip()
                    to_save.append(line.strip())
                if "HESSFREQ" in line and "=" in line:
                    fcfile = line.split("=")[-1].strip()
                    to_save.append(line.strip())
       
        global_save_text.append(to_save)
        return outfile,fcfile


    def load_q2rfc(fcfile):
        global global_natom
        Ry2Hartree = 0.5 
        natom = global_natom

        fm  = []
        flag = "   1   1   1   1"
        for i in range(3*natom):
            fm.append([])
        for i in range(3*natom):
          for j in range(3*natom):
             fm[i].append(0)


        container = []
        lab = 0
        with open(fcfile) as f:
            for line in f:
               if flag in line:
                  lab = 1
                  container.append(line)
                  continue
               if lab == 1 and len(line) > 3:
                  container.append(line)
                  continue
        for i in range(int(len(container)/2)):
            indexline = [int(item) for item in container[2*i].split()]
            a = indexline[0]
            b = indexline[1]
            c = indexline[2]
            d = indexline[3]
            ii = (c-1)*3+ a-1
            jj = (d-1)*3+ b-1
            f = float(container[2*i+1].split()[-1])
            fm[ii][jj] = f*Ry2Hartree
        #print(fm)
        return fm



         


    def load_hessfreq(fcfile):
        fm = []
        fmb = []
        with open(fcfile) as f:
           for line in f:
               if len(line) > 3:
                  [fm.append(float(item)) for item in line.split()]  
        natm3 = int(len(fm)**0.5) 
        fmb = [fm[i:i + natm3] for i in range(0, len(fm), natm3)]
        
        return fmb

    def load_phonopyfc(fcfile,prg):
        fm = []
        evPerHartree = 2.72113838565563E+01
        BohrPerAngstrom = 1.88972613288564E+00
        if prg == "vasp":
           CF = (1/evPerHartree)/( BohrPerAngstrom*BohrPerAngstrom ) #from eV/(Angstrom^2) to Hartree/(Bohr^2) 
        HartreePerRy = 0.5
        if prg == "qe":
           CF = HartreePerRy
        if prg == "cp2k":
           CF = 1/BohrPerAngstrom # not sure


        counter = 0
        with open(fcfile) as f:
           for line in f:
              counter = counter + 1
              if counter == 1:
                 natom = int(line.split()[0]) #int(line)
                 for i in range(3*natom):
                     fm.append([])
                 for i in range(3*natom):
                    for j in range(3*natom):
                        fm[i].append(0)
              if (counter%4) == 2:
                  indi = int(line.split()[0])
                  indj = int(line.split()[1])
              if (counter%4) == 3:
                  xx = float(line.split()[0]) 
                  xy = float(line.split()[1])
                  xz = float(line.split()[2]) 
              if (counter%4) == 0:
                  yx = float(line.split()[0]) 
                  yy = float(line.split()[1])
                  yz = float(line.split()[2]) 
              if (counter%4) == 1 and counter != 1:
                  zx = float(line.split()[0]) 
                  zy = float(line.split()[1])
                  zz = float(line.split()[2])
                  fm[(indi-1)*3+1-1][(indj-1)*3+1-1] = xx *CF
                  fm[(indi-1)*3+1-1][(indj-1)*3+2-1] = xy *CF
                  fm[(indi-1)*3+1-1][(indj-1)*3+3-1] = xz *CF

                  fm[(indi-1)*3+2-1][(indj-1)*3+1-1] = yx *CF
                  fm[(indi-1)*3+2-1][(indj-1)*3+2-1] = yy *CF
                  fm[(indi-1)*3+2-1][(indj-1)*3+3-1] = yz *CF

                  fm[(indi-1)*3+3-1][(indj-1)*3+1-1] = zx *CF
                  fm[(indi-1)*3+3-1][(indj-1)*3+2-1] = zy *CF
                  fm[(indi-1)*3+3-1][(indj-1)*3+3-1] = zz *CF


        return  fm 
   
    def load_qe_xyz(inpfile,obj):
        global global_mass
        global global_save_text  

        to_save_1 = []
        to_save_2 = []
    
        natom = 0 
        with open(inpfile) as f:
            for line in f:
                if "nat" in line:
                    natom = int(line.split("=")[-1])
        f.close()
        if natom == 0:
            print("Error: 'nat' not found") 
            return
        lab = 0
        ct = 0
        ctb = 0
        v=[]
        av = []
        elem = []
        xyz = []
        ifrac = 0
        with open(inpfile) as f:
            for line in f:
                if len(line) > 3:
                    if "CELL_PARAMETERS" in line.upper() and "angstrom" in line.lower():
                       lab = 1
                       continue
                    if lab > 3: 
                       lab =-1
                       #########continue
                    if lab >=1 :
                       a = line.split()
                       av.append( [float(item) for item in a] )
                       astr = a[0]+","+a[1]+","+a[2]
                       v.append(astr)
                       lab = lab + 1
                       continue
                    if "ATOMIC_POSITIONS" in line.upper() and "angstrom" in line.lower():
                       ct = 1
                       continue
                    if ct > natom:
                       ct = -1
                       continue
                    if ct >= 1:
                       at = line.split()[0]
                       global_mass.append(get_mass(at))
                       elem.append(at)
                       ct = ct + 1
                       xyz.append(line)

                    if "ATOMIC_POSITIONS" in line.upper() and "crystal" in line.lower():  
                       ctb = 1
                       ifrac = 1
                       continue
                    if ctb > natom:
                       ctb = -1
                       continue
                    if ctb >= 1:
                       at = line.split()[0]
                       global_mass.append(get_mass(at))
                       elem.append(at)
                       ctb = ctb + 1
                       xyz.append([float(item) for item in line.split()[1:4]])

        if ifrac == 1:
            r=[]
            #print(av)
            #print(xyz)
            for i in range(natom): 
                tm1=list_add(list_scale(av[0],xyz[i][0]),list_scale(av[1],xyz[i][1]))  
                tm2=list_add(tm1,  list_scale(av[2],xyz[i][2]))
                r.append(tm2)


        f1 = open("qe_geom_filename.xyz","w")
        f1.write(str(natom)+"\ntitle\n")
        for i in range(natom):
            if ifrac == 0:
               f1.write(xyz[i])
               to_save_1.append(xyz[i].strip())
            if ifrac == 1:
               f1.write(elem[i]+" "+str(r[i][0])+" "+str(r[i][1])+" "+str(r[i][2])+"\n") 
               to_save_1.append(elem[i]+" "+str(r[i][0])+" "+str(r[i][1])+" "+str(r[i][2]))

        f1.close()
        cmd.load("qe_geom_filename.xyz",obj)
        os.remove("qe_geom_filename.xyz")
        global_save_text.append(to_save_1)


        form.input_v1.setText(v[0])
        form.input_v2.setText(v[1])
        form.input_v3.setText(v[2])
        to_save_2.append(v[0])
        to_save_2.append(v[1])
        to_save_2.append(v[2])
        global_save_text.append(to_save_2)


        return natom,elem 

    def load_cp2k_cell(inpfile):
        global global_save_text

        to_save = []

        flag = "&CELL"
        lab = 0
        iread = 0
        v=[]
        with open(inpfile) as f:
            for line in f:
                if len(line) >3 and flag.lower() in line.lower():
                    if iread == 0:
                       lab = 1 
                    continue
                if lab > 3:
                   lab = -1
                   continue
                if lab >=1:
                   iread = 1 
                   v.append([float(item) for item in line.split()[1:4]])
                   lab = lab + 1
                   continue
        v1s = str(v[0][0])+","+str(v[0][1])+","+str(v[0][2])
        v2s = str(v[1][0])+","+str(v[1][1])+","+str(v[1][2])
        v3s = str(v[2][0])+","+str(v[2][1])+","+str(v[2][2])
        form.input_v1.setText(v1s)
        form.input_v2.setText(v2s)
        form.input_v3.setText(v3s)
        to_save.append( v1s )
        to_save.append( v2s )
        to_save.append( v3s )
        global_save_text.append(to_save)

    
        return

    def load_fchk(fchkfile,obj):
        import numpy as np
        global global_mass
        global global_save_text

        b2a = 0.529177
        elem = []
        natom = 0
        xyzl = []
        hessl = []

        lab1 = 0
        lab2 = 0
        lab3 = 0
        with open(fchkfile) as f:
            for line in f:
                if "Number of atoms" in line:
                   natom = int(line.split()[-1])
                   nrowz = int(natom / 6) + 1
                   if natom % 6 == 0:
                      nrowz = nrowz - 1
                   nrowc = int(3*natom / 5) + 1
                   if (3*natom) % 5 == 0:
                      nrowc = nrowc - 1 
                   continue
                if "Atomic numbers" in line:
                   lab1 = 1
                   continue
                if lab1 >= 1 and lab1 <= nrowz:
                   lab1 = lab1 + 1
                   [elem.append(get_symbol(int(item))) for item in line.split()]
                   [global_mass.append(get_mass(get_symbol(int(item)))) for item in line.split()]   
                   continue
                if "Current cartesian coordinates" in line:
                   lab2 = 1
                   continue
                if lab2 >= 1 and lab2 <= nrowc:
                   lab2 = lab2 + 1
                   [xyzl.append(b2a*float(item)) for item in line.split()]
                   continue
                if "Cartesian Force Constants" in line:
                   lab3 = 1
                   nrowh = int(int(line.split()[-1]) / 5) + 1
                   if int(line.split()[-1]) % 5 == 0:
                      nrowh = nrowh - 1 
                   continue
                if lab3 >= 1 and lab3 <= nrowh:
                   lab3 = lab3 + 1 
                   [hessl.append(float(item)) for item in line.split()]
                   continue

                #
        #
        to_save_1 = [] 
        f1 = open("fchk_geom_filename.xyz","w")
        f1.write(str(len(elem))+"\ntitle\n")
        for i in range(len(elem)):
            f1.write(elem[i]+" "+ str(xyzl[3*i])+" "+ str(xyzl[3*i+1])+" "+ str(xyzl[3*i+2])+"\n")
            to_save_1.append(elem[i]+" "+ str(xyzl[3*i])+" "+ str(xyzl[3*i+1])+" "+ str(xyzl[3*i+2]))
        f1.close()
        cmd.load("fchk_geom_filename.xyz",obj)
        os.remove("fchk_geom_filename.xyz")
        global_save_text.append(to_save_1) # xyz
        global_save_text.append([]) # cell
               
        #
        size_X = 3*natom
        X = np.zeros((size_X,size_X))
        X[np.tril_indices(X.shape[0], k = 0)] = np.array(hessl)
        X = X + X.T - np.diag(np.diag(X))
        X = X.tolist()

        return natom,elem,X

        

    def load_unimovib(umvfile,obj):
        global global_mass
        global global_save_text
       
        to_save_1 = []

        b2a = 0.529177
        elem = []
        natom = 0
        mass = []
        xyzl = []
        fm = []

        lab = 0
        #print(umvfile)
        with open(umvfile) as f:
            for line in f :
              if "NATM" in line:
                 lab = 1
                 continue
              if lab == 1:
                 natom = int(line) 
                 lab = 2
                 continue
        f.close()
      
        nrowz = int(natom / 5) + 1
        if natom % 5 == 0:
           nrowz = nrowz - 1
        lab = 0   
        with open(umvfile) as f:
            for line in f :
              if "ZA" in line:
                 lab = 1
                 continue
              if lab >= 1 and lab <= nrowz:
                 #natom = int(line)
                 
                 z = [int(eval(item.replace("D","E"))) for item in line.split()]
                 [elem.append(get_symbol(item)) for item in z] 
                 [global_mass.append(get_mass(get_symbol(item))) for item in z]

                 lab = lab+1
                 continue
        f.close()
             
        nrow = int((3*natom) /5)+1
        if (3*natom) % 5 == 0:
           nrow = nrow - 1
        lab = 0
        with open(umvfile) as f:
            for line in f:
              if "XYZ" in line:
                 lab = 1
                 continue
              if lab >= 1 and lab <= nrow:
                 x = [b2a*float(eval(item.replace("D","E"))) for item in line.split()]
                 [xyzl.append(item) for item in x]
                 lab = lab + 1
                 continue
        f.close()


        f1 = open("unimovib_geom_filename.xyz","w")
        f1.write(str(len(elem))+"\ntitle\n")
        for i in range(len(elem)):
            f1.write(elem[i]+" "+ str(xyzl[3*i])+" "+ str(xyzl[3*i+1])+" "+ str(xyzl[3*i+2])+"\n")
            to_save_1.append(elem[i]+" "+ str(xyzl[3*i])+" "+ str(xyzl[3*i+1])+" "+ str(xyzl[3*i+2]))
        f1.close()
        cmd.load("unimovib_geom_filename.xyz",obj)
        os.remove("unimovib_geom_filename.xyz")
        global_save_text.append(to_save_1) # xyz
        global_save_text.append([]) # cell
        

        select_mol() # radio

             
        nrow = int((3*natom*3*natom) /5)+1
        if (3*natom*3*natom) % 5 == 0:
           nrow = nrow - 1
        lab = 0
        hessl = []
        with open(umvfile) as f:
            for line in f:
              if "FFX" in line:
                 lab = 1
                 continue
              if lab >= 1 and lab <= nrow:
                 x = [float(eval(item.replace("D","E"))) for item in line.split()]
                 [hessl.append(item) for item in x]
                 lab = lab + 1
                 continue
        f.close()
        
        for i in range(3*natom):
            fm.append([])    
        for i in range(3*natom):
            for j in range(3*natom):
                fm[i].append(0)
        for i in range(3*natom):
            for j in range(3*natom):
                fm[i][j]  =   hessl[3*natom*i + j]     

        
        return natom,elem,fm
      

    def load_castep_xyz(outfile,obj):
        global global_mass
        global global_save_text

        elem = []
        natom = 0

        flag_1 = "Real Lattice"
        lab = 0
        v=[]
        with open(outfile) as f:
            for line in f :
              if "Total number of ions in cell =" in line:
                 natom = int(line.split()[-1])
                 continue
              if flag_1 in line:
                 lab = 1
                 continue
              if lab >= 1 and lab <= 3:
                 v.append([float(item) for item in line.split()[0:3]])
                 lab = lab + 1
                 continue
        f.close()

        
        lab = 0
        flag_2 = "Element    Atom"
        xyz_frac = []
        with open(outfile) as f:
            for line in f:
               if flag_2 in line: 
                  lab = 1
                  continue
               if lab >= 1 and lab < 3 :
                  lab = lab + 1
                  continue
               if lab >= 3 and lab < (3+natom) and len(global_mass) < natom :
                  #print(line) 
                  elem.append(line.split()[1])
                  e = line.split()[1]
                  global_mass.append(get_mass(e))
                  xyz_frac.append([float(item) for item in line.split()[3:6]])

                  lab = lab + 1
                  continue
        f.close()

        xyz_car = []
        for i in range(natom):
               tm1=list_add(list_scale(v[0],xyz_frac[i][0]),list_scale(v[1],xyz_frac[i][1]))
               r=list_add(tm1,list_scale(v[2],xyz_frac[i][2]))
               xyz_car.append(r)

        to_save_1 = []
        f1 = open("castep_geom_filename.xyz","w")
        f1.write(str(natom)+"\ntitle\n")
        for i in range(natom):
            f1.write(elem[i]+" "+str(xyz_car[i][0])+" "+str(xyz_car[i][1])+" "+str(xyz_car[i][2])+"\n")
            to_save_1.append(elem[i]+" "+str(xyz_car[i][0])+" "+str(xyz_car[i][1])+" "+str(xyz_car[i][2]))
        f1.close()
        cmd.load("castep_geom_filename.xyz",obj)
        os.remove("castep_geom_filename.xyz")
        global_save_text.append(to_save_1)

        v1s = str(v[0][0])+","+str(v[0][1])+","+str(v[0][2])
        v2s = str(v[1][0])+","+str(v[1][1])+","+str(v[1][2])
        v3s = str(v[2][0])+","+str(v[2][1])+","+str(v[2][2])
        form.input_v1.setText(v1s)
        form.input_v2.setText(v2s)
        form.input_v3.setText(v3s)
        to_save = []
        to_save.append( v1s )
        to_save.append( v2s )
        to_save.append( v3s )
        global_save_text.append(to_save)


        # read in dynamical matrix and obtain hessian 
        flag = "Dynamical matrix"
        nhit = 0
        with open(outfile) as f:
            for line in f:
                if flag in line:
                   nhit = nhit + 1 
        f.close()
        ihit = 0
        lab = 0 
        nrow = int((3*natom) / 6) + 1 
        if (3*natom) % 6 == 0:
           nrow = nrow -1 

        dm = []
        dml = []
        for i in range(3*natom):
            dm.append([])
        for i in range(3*natom):
            for j in range(3*natom):
                dm[i].append(0)



        with open(outfile) as f:
            for line in f:
                if flag in line:
                   ihit = ihit + 1 
                   continue
                if ihit == nhit:
                   lab = lab + 1
                   ihit = -1
                   continue
                if lab >= 1 and lab <= (3*natom*nrow):
                   ss = line.split()
                   for s in ss:
                       if "." in s:
                           dml.append(float(s))

                   lab  = lab + 1
                   continue
        f.close()
        for i in range(3*natom):
            for j in range(3*natom):
                dm[i][j] = dml[i*3*natom+j]
        #print(dm)        
          
        for i in range(3*natom):
            for j in range(3*natom):
                ni = int(i / 3)
                nj = int(j / 3)
                dm[i][j] = dm[i][j]/(au2wn*au2wn)*(get_mass(elem[ni])*get_mass(elem[nj]))**0.5 




        return natom,elem,dm

    def load_cp2k_nexyz(inpfile,obj):
        global global_mass
        global global_save_text
       
        to_save_1 = []

        elem = []
        natom = 0
        mass = []
        flag_1 = "&COORD"
        lab = 0
        xyz_str = []
        with open(inpfile) as f:
             for line in f:
               if flag_1 in line.upper():
                  lab = 1
                  continue
               if lab == 1:
                  if "END" in line.upper():
                     lab = 2
                     break
                  else:
                     natom = natom + 1  
                     e = line.split()[0]
                     global_mass.append(get_mass(e))
                     elem.append(e)
                     xyz_str.append(line.strip())
                     continue
        f.close()

        f1 = open("cp2k_geom_filename.xyz","w")
        f1.write(str(len(elem))+"\ntitle\n")
        for i in range(len(elem)):
            f1.write(xyz_str[i]+"\n")
            to_save_1.append(xyz_str[i])
        f1.close()
        cmd.load("cp2k_geom_filename.xyz",obj)
        os.remove("cp2k_geom_filename.xyz")
        global_save_text.append(to_save_1)

        
        return natom,elem


    def load_cp2k_xyz(outfile,obj):
        global global_mass
        global global_save_text

        to_save_1 = []

        flag_1 = "Atom  Kind  Element       X"
        elem = []
        natom = 0
        lab = 0
        istart=0
        coor_str=[]
        mass = []
        with open(outfile) as f:
            for line in f :
              if flag_1 in line:
                 lab = 1
                 continue
              #if lab == 1:
              #   if "1" not in line: 
              #      lab = 2
              #   else:
              #
              #   continue
              if lab == 1:
                 if ("1" not in line) and istart==0 :
                    continue

                 if len(line.strip()) != 0:
                    istart=1 
                    natom = natom + 1 
                    e = line.split()[2]
                    global_mass.append(get_mass(e))
                    elem.append(e)
                    x = str(float(line.split()[4]))
                    y = str(float(line.split()[5]))
                    z = str(float(line.split()[6]))
                    coor_str.append([x,y,z])
                    m = float(line.split()[8])
                    for i in range(3): 
                       mass.append(m)
                    continue
                 else:
                    if istart==1: 
                       lab = 3
                       break
        f.close()

        f1 = open("cp2k_geom_filename.xyz","w")
        f1.write(str(len(elem))+"\ntitle\n")
        for i in range(len(elem)):
            f1.write(elem[i]+" "+coor_str[i][0]+" "+coor_str[i][1]+" "+coor_str[i][2]+"\n")
            to_save_1.append(elem[i]+" "+coor_str[i][0]+" "+coor_str[i][1]+" "+coor_str[i][2])
        f1.close()
        cmd.load("cp2k_geom_filename.xyz",obj)
        os.remove("cp2k_geom_filename.xyz")
        global_save_text.append(to_save_1)

 
        with open(outfile) as f:
            for line in f:
                if "CELL| Vector a [angstrom]" in line:
                   v = line.split()[4]+","+line.split()[5]+","+line.split()[6]  
                   #form.input_v1.setText(v)
                if "CELL| Vector b [angstrom]" in line:
                   v = line.split()[4]+","+line.split()[5]+","+line.split()[6]  
                   #form.input_v2.setText(v)
                if "CELL| Vector c [angstrom]" in line:
                   v = line.split()[4]+","+line.split()[5]+","+line.split()[6]  
                   #form.input_v3.setText(v)
                  


        f.close()
        lab = 0
        container = []
        flag = "VIB| Hessian in cartesian coordinates"
        with open(outfile) as f:
            for line in f:
                if flag in line:
                   lab = 1
                   continue
                if "VIB| Cartesian" in line: #or "VIB| Frequencies after" in line:
                   lab = -1

                if lab == 1:
                   container.append(line)
                   continue

        f.close()
        #print(container[0])
        #print(container[1])
        #print(container[2])

        #print(container[-3])
        #print(container[-2])
        #print(container[-1])
        #print("natom:", natom)
   
        fm = []
        fc = 1.660538782E-27/9.10938215E-31*1E-6 
        for i in range(3*natom):
            fm.append([])
        for i in range(3*natom):
            for j in range(3*natom):
                   fm[i].append(0)
        #
        for i in range(3*natom):
            for j in range(3*natom):
                
                block_row = int(i/5)
                block_col = int(i%5)

                n = 3*natom + 2
                row = block_row * n + 2
                col = block_col + 2 
             
                fm[j][i] = float(container[row+j].split()[col])
                fm[j][i] = fm[j][i]*(mass[i]*mass[j])**0.5 *fc

        #print(fm)

        #
        return natom,elem,fm


    def load_crystal17_xyz(outfile,obj):
        global global_mass
        global global_save_text

        to_save_1 = []
        to_save_2 = []

        imol = 0
        # molecular calculation
        flag_0 = "MOLECULAR CALCULATION"
        flag_1m= "ATOM AT. N.              COORDINATES"

        with open(outfile) as f:
            for line in f:
              if flag_0 in line:
                 imol = 1
                 break
        f.close()



        # non-molecular calculation
        flag_1 = "      ATOM          X(ANGSTROM)         Y(ANGSTROM)         Z(ANGSTROM)"
        flag_2 = "DIRECT LATTICE VECTORS CARTESIAN COMPONENTS (ANGSTROM)"

        elem = []
        lab = 0
        coor_str = []
        natom = 0
        with open(outfile) as f:
            for line in f:
              if flag_1 in line:
                 lab = 1
                 continue
              if lab == 1:
                 lab = 2
                 continue
              if lab == 2:
                 if len(line.strip()) != 0:
                    natom = natom + 1 
                    e = line.split()[2]
                    elem.append(e)
                    global_mass.append(get_mass(e))
                    x = str(float(line.split()[3]))
                    y = str(float(line.split()[4]))
                    z = str(float(line.split()[5]))
                    coor_str.append([x,y,z])
                    continue
                 else:
                    lab = 3
                    break
        f.close()

        # mol
        if imol == 1:
           with open(outfile) as f:
               for line in f:
                  if flag_1m in line:
                     lab = 1
                     continue
                  if lab == 1:
                     if len(line.strip()) != 0:
                        natom = natom + 1
                        iz = int(line.split()[1])
                        e = get_symbol(iz)
                        elem.append(e)
                        global_mass.append(get_mass(e))
                        x = str(float(line.split()[2]))
                        y = str(float(line.split()[3]))
                        z = str(float(line.split()[4]))
                        coor_str.append([x,y,z])
                        continue
                     else:
                        lab = 2
                        break
           f.close()   
 


        f1 = open("crystal17_geom_filename.xyz","w")
        f1.write(str(len(elem))+"\ntitle\n")
        for i in range(len(elem)):
            f1.write(elem[i]+" "+coor_str[i][0]+" "+coor_str[i][1]+" "+coor_str[i][2]+"\n")
            to_save_1.append(elem[i]+" "+coor_str[i][0]+" "+coor_str[i][1]+" "+coor_str[i][2])

        f1.close()
        cmd.load("crystal17_geom_filename.xyz",obj)
        os.remove("crystal17_geom_filename.xyz")
        global_save_text.append(to_save_1)
        
        lab = 0
        v=[]
        with open(outfile) as f:
            for line in f:
              if flag_2 in line:
                 lab = 1
                 continue
              if lab == 1:
                 lab = 2
                 continue
              if lab == 2:
                 if len(line.strip()) != 0:
                    x = str(float(line.split()[0]))
                    y = str(float(line.split()[1]))
                    z = str(float(line.split()[2]))
                    v.append(x+","+y+","+z)
                    continue
                 else:
                    lab = 3
                    break
        f.close()
        
        if imol == 1:
           v=["500.0,0.0,0.0","0.0,500.0,0.0","0.0,0.0,500.0"] 

        form.input_v1.setText(v[0])
        form.input_v2.setText(v[1])
        form.input_v3.setText(v[2])
        to_save_2.append(v[0])
        to_save_2.append(v[1])
        to_save_2.append(v[2])
        global_save_text.append(to_save_2)

        
        #
        return natom,elem


    def load_poscar(posfile,obj):
        global global_mass 
        global global_save_text
        
        to_save_1 = []
        to_save_2 = []

        f = 1.0
        dat=[]
        with open(posfile) as f:
            for i,line in enumerate(f):
                if i == 1:
                   f = float(line.strip()) 
                if i == 2:
                   v1 = line.split()
                   v1 = [f*float(a) for a in v1]
                   v1s = str(v1[0])+","+str(v1[1])+","+str(v1[2])
                   form.input_v1.setText(v1s)
                   to_save_2.append(v1s)
                if i == 3:
                   v2 = line.split()
                   v2 = [f*float(a) for a in v2]
                   v2s = str(v2[0])+","+str(v2[1])+","+str(v2[2])
                   form.input_v2.setText(v2s)
                   to_save_2.append(v2s)
                if i == 4:
                   v3 = line.split()
                   v3 = [f*float(a) for a in v3]
                   v3s = str(v3[0])+","+str(v3[1])+","+str(v3[2])
                   form.input_v3.setText(v3s)
                   to_save_2.append(v3s)
                if i == 5:
                   elem = line.split()  
                if i == 6:
                   count = line.split()
                   count = [int(a) for a in count]
                   natom = sum(count) 

                if i == 7:
                   if line[0] == "C" or line[0] == "c":
                      t= 1

                   if line[0] == "D" or line[0] == "d":    
                      t= 2
                if i>=8 and i<(8+natom):
                   d = line.split()
                   d = [float(a) for a in d]
                   dat.append(d)
                   #print(d)
        #
        xyz = []
        if t == 1:
           for i in range(natom):
               dat[i] = [f*b for b in dat[i]]      
           xyz = dat
        if t == 2:
           for i in range(natom):
               tm1=list_add(list_scale(v1,dat[i][0]),list_scale(v2,dat[i][1]))
               r=list_add(tm1,list_scale(v3,dat[i][2]))
               xyz.append(r)
        
        atoms = []
        for i in range(len(elem)):
            c = count[i]
            for j in range(c):
                atoms.append(elem[i])
                atm = get_mass(elem[i])
                global_mass.append(atm)

        
        f1 = open("vasp_geom_filename.xyz","w")
        f1.write(str(natom)+"\ntitle\n")
        for i in range(natom):
            f1.write(atoms[i]+" "+str(xyz[i][0])+" "+str(xyz[i][1])+" "+str(xyz[i][2])+"\n")
            to_save_1.append(atoms[i]+" "+str(xyz[i][0])+" "+str(xyz[i][1])+" "+str(xyz[i][2]))
        f1.close()
        cmd.load("vasp_geom_filename.xyz",obj)
        os.remove("vasp_geom_filename.xyz")
        global_save_text.append(to_save_1)
        global_save_text.append(to_save_2)


        return natom,atoms



    def loadxyz():
        global global_natom
        global global_elem
        global global_fm
        global global_minv
        global global_save_text

        xyz_file_path = form.input_geomname.text().strip()
        if len(xyz_file_path) == 0:
           return
        prg = form.list_program.currentIndex()  
        pymol_init()

        if prg == 6: # UniMoVib
           global_save_text.append("UniMoVib generic interface")
           global_save_text.append(xyz_file_path)

           umvf = read_unimovib(xyz_file_path)
           dirp = os.path.dirname(xyz_file_path)
           umvf = os.path.join(dirp,umvf)

           if (not os.path.exists(umvf)):
              QtWidgets.QMessageBox.critical(None, 'Error', "Something goes wrong with "\
              +xyz_file_path+\
              "  Please check the path for UMV data file!")
              return 
           global_natom,global_elem,global_fm = load_unimovib(umvf,"geom")
        
        if prg == 5: # gaussian/qchem
           global_save_text.append("Gaussian/Q-Chem")
           global_save_text.append(xyz_file_path)
         

           global_save_text.append(["please see the formatted checkpoint file."])
             
           # read in formatted checkpoint file directly?
           if (not os.path.exists(xyz_file_path)):
              QtWidgets.QMessageBox.critical(None, 'Error', "Something goes wrong with "\
              +xyz_file_path+\
              "  Please check the path for Gaussian/Q-Chem formatted checkpoint file!")
              return
           global_natom,global_elem,global_fm = load_fchk(xyz_file_path,"geom")

           #return 

        if prg == 4: # castep
           global_save_text.append("CASTEP")
           global_save_text.append(xyz_file_path)
           
           outf = read_castep(xyz_file_path)
           dirp = os.path.dirname(xyz_file_path)
           #clf = os.path.join(dirp,clf)
           outf = os.path.join(dirp,outf)
           #print(outf)
           if (not os.path.exists(outf))  :
              QtWidgets.QMessageBox.critical(None, 'Error', "Something goes wrong with "\
              +xyz_file_path+\
              "  Please check the path for CASTEP output file!")
              return
           global_natom,global_elem,global_fm = load_castep_xyz(outf,"geom")
           
           #return

        if prg == 3: # cp2k
           global_save_text.append("CP2K")
           global_save_text.append(xyz_file_path)

           mode,fcf,inpf = read_cp2k(xyz_file_path)
           dirp = os.path.dirname(xyz_file_path)
           fcf  = os.path.join(dirp,fcf)
           inpf = os.path.join(dirp,inpf)
           if (not os.path.exists(fcf)) or (not os.path.exists(inpf))  :
              QtWidgets.QMessageBox.critical(None, 'Error', "Something goes wrong with "\
              +xyz_file_path+\
              "  Please check the path for CP2K INPUT and OUTPUT/FCPHONOPY file!")
              return
           if mode == "native":
              global_natom,global_elem,global_fm = load_cp2k_xyz(fcf,"geom")
           if mode == "phonopy":
              global_natom,global_elem = load_cp2k_nexyz(inpf,"geom")   
              global_fm = load_phonopyfc(fcf,"cp2k")
              pass 
           load_cp2k_cell(inpf)

        
        if prg == 2: # QE
           global_save_text.append("Quantum ESPRESSO") 
           global_save_text.append(xyz_file_path)
           mode,inpf,fcf = read_qe(xyz_file_path)
           dirp = os.path.dirname(xyz_file_path)
           inpf = os.path.join(dirp,inpf)
           fcf  = os.path.join(dirp,fcf)
           if (not os.path.exists(inpf)) or (not os.path.exists(fcf)):
              QtWidgets.QMessageBox.critical(None, 'Error', "Something goes wrong with "\
              +xyz_file_path+\
              "  Please check the path for QE input geometry file and q2r force constant file!")       
              return 
 
           if mode == 'dfpt':
              global_natom,global_elem = load_qe_xyz(inpf,"geom")
              global_fm = load_q2rfc(fcf)
           if mode == 'phonopy':
              global_natom,global_elem = load_qe_xyz(inpf,"geom")
              global_fm = load_phonopyfc(fcf,"qe")
              pass 


        if prg == 1: # crystal17
           global_save_text.append("CRYSTAL") 
           global_save_text.append(xyz_file_path)
           outf,fcf = read_crystal(xyz_file_path)
           dirp = os.path.dirname(xyz_file_path)
           outf = os.path.join(dirp,outf)
           fcf  = os.path.join(dirp,fcf)
           if (not os.path.exists(outf)) or (not os.path.exists(fcf)):
              QtWidgets.QMessageBox.critical(None, 'Error', "Something goes wrong with "\
              +xyz_file_path+\
              "  Please check the path for CRYSTAL OUTPUT file and HESSFREQ file!")       
              return 
 
           global_natom,global_elem = load_crystal17_xyz(outf,"geom")  
           global_fm = load_hessfreq(fcf)

        if prg == 0: # vasp
           global_save_text.append("VASP") 
           global_save_text.append(xyz_file_path)
           posf,fcf = read_vasp(xyz_file_path)
           dirp = os.path.dirname(xyz_file_path)
           posf = os.path.join(dirp,posf)
           fcf  = os.path.join(dirp,fcf)
           
           if (not os.path.exists(posf)) or (not os.path.exists(fcf)):
              QtWidgets.QMessageBox.critical(None, 'Error', "Something goes wrong with "\
              +xyz_file_path+\
              "  Please check the path for POSCAR/CONCAR file and FORCE_CONSTANT file!")       
              return 
           #print(posf,fcf)
           
           global_natom,global_elem = load_poscar(posf,"geom")
           global_fm = load_phonopyfc(fcf,"vasp")

        #print("mass:",global_mass)

        global_minv = calc_minv(global_mass)

 
        set_valence_obj("geom",1)
        set_style("geom",1,"0.13","0.05")
        form.tableWidget.setRowCount(0)
        form.button_remove.setEnabled(False) 


    def loadxyz_temp():
        global global_natom
        global global_elem 
        xyz_file_path = form.input_geomname.text().strip()
        if len(xyz_file_path) == 0:
           return
        pymol_init()

        cmd.load(xyz_file_path,"geom")
        xyz_coor = cmd.get_model('geom', 1).get_coord_list()
        global_natom = len(xyz_coor)

        model_o = cmd.get_model('geom', 1)
        #elem = []
        for at in model_o.atom:
            global_elem.append(at.name)



        set_valence_obj("geom",1)

        set_style("geom",1,"0.13","0.05")
        form.tableWidget.setRowCount(0)#clear up the table 



   

    def pymol_init():
        cmd.bg_color("white")
        cmd.set("ray_shadow","0") 
        cmd.set("orthoscopic") # don't use perspective
        cmd.set("auto_zoom","off")
        cmd.delete("geom")

        cmd.set("antialias_shader","2")
        cmd.set("line_smooth","1")
        cmd.set("depth_cue","1")
        cmd.set("specular","1.00000")
        cmd.set("surface_quality","1")
        cmd.set("cartoon_sampling","14")
        cmd.set("ribbon_sampling","10")
        cmd.set("transparency_mode","2")
        cmd.set("use_shaders","1")
        cmd.set("cartoon_use_shader","1")
        cmd.set("cgo_use_shader","1")
        cmd.set("dash_use_shader","1")
        cmd.set("dot_use_shader","1")
        cmd.set("line_use_shader","1")
        cmd.set("mesh_use_shader","1")
        cmd.set("nb_spheres_use_shader","1")
        cmd.set("nonbonded_use_shader","1")
        cmd.set("ribbon_use_shader","1")
        cmd.set("sphere_use_shader","1")
        cmd.set("stick_use_shader","1")
        cmd.set("surface_use_shader","1")
        cmd.set("render_as_cylinders","1")
        cmd.set("alignment_as_cylinders","1")
        cmd.set("cartoon_nucleic_acid_as_cylinders","1")
        cmd.set("dash_as_cylinders","1")
        cmd.set("line_as_cylinders","1")
        cmd.set("mesh_as_cylinders","1")
        cmd.set("nonbonded_as_cylinders","1")
        cmd.set("ribbon_as_cylinders","1")
        cmd.set("stick_as_cylinders","1")
        cmd.set("dot_as_spheres","1")
        cmd.set("stick_ball","0")
        cmd.set("sphere_mode","9")
        cmd.set("nb_spheres_quality","3")

        return 

        
    def set_style(objname,n,sp_scale,ln_rad):
        
        if n == 1:
           cmd.show_as("lines",objname)
           cmd.show("spheres",objname)
           cmd.set("sphere_scale",sp_scale,objname)
           set_atomcolor()
           cmd.set("line_radius", ln_rad,objname) 

        else:
           return 

    def set_atomcolor():
        cmd.color("grey6","elem C")
        cmd.color("red", "elem O")
        cmd.color("white", "elem H")

        cmd.color("aluminum", "elem Al")
        cmd.color("bismuth", "elem Bi")
        cmd.color("boron","elem B")
        cmd.color("bromine", "elem Br")
        cmd.color("potassium","elem K")
        cmd.color("rhenium","elem Re")


    def set_valence_obj(obj_name,insertQ): 
        model_o = cmd.get_model(obj_name, 1)
        elem = []
        for at in model_o.atom:
            elem.append(at.name)
        xyz_list = cmd.get_model(obj_name, 1).get_coord_list()
        set_valence(elem,xyz_list,obj_name,insertQ)

    #def calc_dis(l1,l2):
    #    d = 0.0
    #    for i in range(3):
    #       d = d + (l1[i]-l2[i])**2
    #    d = d**0.5
    #    return d 


    def judge_valence(lis2,dis):
        val = 0.0
        if lis2 == ["C","C"]:
           val = 0.0
           if dis < 1.64:
              val = 1.0
           if dis < 1.45:
              val = 1.5
           if dis < 1.39:
              val = 2.0
           if dis < 1.25:
              val = 3.0

        if lis2 == ["C","O"] or lis2 == ["O","C"]:
           val = 0.0
           if dis < 1.52:
              val = 1.0
           if dis < 1.35: 
              val = 1.5
           if dis < 1.29:
              val = 2.0
           if dis < 1.16:
              val = 3.0
        
        if lis2 == ["C","N"] or lis2 == ["N","C"]:
           val = 0.0
           if dis < 1.56:
              val = 1.0
           if dis < 1.39:
              val = 1.5 
           if dis < 1.33:
              val = 2.0
           if dis < 1.20:
              val = 3.0

        if lis2 == ["H","H"]:
           val = 0.0
           if dis <= 1.0:
              val = 1.0

        # user can add more valence information when necessary.
        if lis2 == ["Na","Cl"] or lis2 == ["Cl","Na"]:
           val = 0.0
           if dis <= 2.50:
              val = 1.0  

        if lis2 == ["Re","H"] or lis2 == ["H","Re"]:
           val = 0.0
           if dis <= 1.72:
              val = 1.0  


        return val    



    def set_valence(elem,xyz_list,obj_name,insertQ):
        global global_delocalized_bonds_list

        insert_flag = 0
        if insertQ == 1:
           insert_flag = 1 

        # this function is to set up correct valence  
        for i in range(len(elem)):
           for j in range(len(elem)):
               if j > i:
                   a = elem[i]
                   b = elem[j]
                   a_xyz = xyz_list[i]
                   b_xyz = xyz_list[j] 
                   d = calc_dis(a_xyz,b_xyz)
                   #print("flag1",[a,b],d)
                   valence = judge_valence([a,b],d)
                   #print("flag2")
                   if valence != 0.0: # and valence != 1.0:
                      #print(a,b)
                      cmd.delete("pka")
                      cmd.delete("pkb")
                      pkaa="id "+str(i+1)+" and "+" model "+obj_name
                      pkbb="id "+str(j+1)+" and "+" model "+obj_name
                      cmd.select("pka","id "+str(i+1)+" and "+" model "+obj_name) # select atom here...
                      cmd.select("pkb","id "+str(j+1)+" and "+" model "+obj_name)
                      if valence == 1.0:
                         #cmd.bond("pka","pkb") # just in case no bond is shown
                         cmd.bond("id "+str(i+1)+" and "+" model "+obj_name,"id "+str(j+1)+" and "+" model "+obj_name)
                         # above line is to correct a bug in cmd.bond 
                         #cmd.valence("1","pka","pkb")
                         # above line is causing another issue
                         cmd.valence("1","id "+str(i+1)+" and "+" model "+obj_name,"id "+str(j+1)+" and "+" model "+obj_name)


                         cmd.set("valence_mode","0")
                         cmd.set("valence_size","0.06")
                      if valence == 2.0:
                         #cmd.valence("2","pka","pkb")
                         cmd.valence("2",pkaa,pkbb)
                         cmd.set("valence_mode","0")
                         cmd.set("valence_size","0.06")
                      if valence == 3.0:
                         cmd.valence("3",pkaa,pkbb) 
                         #cmd.valence("3","pka","pkb")
                         cmd.set("valence_mode","0")
                         cmd.set("valence_size","0.06")
                      if valence == 1.5:
                         cmd.valence("4",pkaa,pkbb) 
                         #cmd.valence("4","pka","pkb")
                         cmd.set("valence_mode","0")
                         cmd.set("valence_size","0.06")
                         if insert_flag == 1:
                            global_delocalized_bonds_list.append([i,j]) # store delocalized bonds 

                      cmd.delete("pka")
                      cmd.delete("pkb") 


    def update_dimension():
        global global_dim
        global global_save_text
        dimension = form.input_dimension.currentText()
        global_dim = dimension


        if dimension == "0":
           form.input_v1.setDisabled(True)
           form.input_v2.setDisabled(True)
           form.input_v3.setDisabled(True)

        if dimension == "1":
           form.input_v1.setEnabled(True)
           form.input_v2.setDisabled(True)
           form.input_v3.setDisabled(True)

        if dimension == "2":
           form.input_v1.setEnabled(True)
           form.input_v2.setEnabled(True)
           form.input_v3.setDisabled(True)
          
        if dimension == "3":
           form.input_v1.setEnabled(True)
           form.input_v2.setEnabled(True)
           form.input_v3.setEnabled(True)

        if len(global_save_text) < 6:
           if len(form.input_geomname.text()) != 0: 
              global_save_text.append( int(dimension) )
        else:
           global_save_text[5] = int(dimension)

          
        return 
      
    def select_mol():
        
        form.input_dimension.setCurrentIndex(0)
        update_dimension()

    def get_min_corner(xyz): # as the origin point of basis vectors 
        x = xyz[0][0]
        y = xyz[0][1]
        z = xyz[0][2]
        for i in range(len(xyz)):
            if x > xyz[i][0]:
               x = xyz[i][0]
            if y > xyz[i][1]:
               y = xyz[i][1]
            if z > xyz[i][2]:
               z = xyz[i][2]
        return [x,y,z] 

    def get_max_corner(xyz): # as the origin point of basis vectors 
        x = xyz[0][0]
        y = xyz[0][1]
        z = xyz[0][2]
        for i in range(len(xyz)):
            if x < xyz[i][0]:
               x = xyz[i][0]
            if y < xyz[i][1]:
               y = xyz[i][1]
            if z < xyz[i][2]:
               z = xyz[i][2]
        return [x,y,z] 




    #def list_minus(a,b): # a-b
    #    c = []
    #    if len(a) != len(b):
    #       return c
    #    for i in range(len(a)):
    #        c.append(a[i]-b[i])
    #    return c


    #def list_add(a,b):
    #    c = []
    #    if len(a) != len(b):
    #       return c
    #    for i in range(len(a)):
    #        c.append(a[i]+b[i])
    #    return c

    def basis_style(color,name):
        cmd.hide("label",name)
        cmd.set("dash_color",color,name)
        cmd.set("dash_gap","0.01",name)
        cmd.set("dash_radius","0.020",name)
        cmd.set("dash_transparency","0.75",name)



    def append_xyz(old,add):
        for i in range(len(add)):
            old.append(add[i])

        return old


    def translate_geom(xyz, vector, direction):
        new_xyz = []
        if direction == "+":
           for i in range(len(xyz)):
               x1 = xyz[i][0] + vector[0]
               y1 = xyz[i][1] + vector[1]
               z1 = xyz[i][2] + vector[2]
               new_xyz.append([x1,y1,z1])
        else:
           for i in range(len(xyz)): 
               x1 = xyz[i][0] - vector[0]
               y1 = xyz[i][1] - vector[1]
               z1 = xyz[i][2] - vector[2]
               new_xyz.append([x1,y1,z1])

        return new_xyz

    def write_xyz(xyz_filename,supercell_xyz,elem):
        try:
           os.remove(xyz_filename)
        except OSError:
           pass
        f1=open(xyz_filename,"w")
        f1.write(str(len(supercell_xyz))+"\n")
        f1.write("title\n")
        for i in range(len(supercell_xyz)):
            j = i%len(elem)
            f1.write( elem[j] + " ")
            f1.write( str(supercell_xyz[i][0]) + " " )
            f1.write( str(supercell_xyz[i][1]) + " " )
            f1.write( str(supercell_xyz[i][2]) + "\n" )
        f1.close()

    #def list_dotp(a,b):
    #    n = len(a)
    #    s =  0.0
    #    for i in range(n):
    #        s = s + a[i]*b[i]
    #    return s 
    
    def hide_outer(n,natom,i,objname):
        cmd.select("pp","id "+str(natom*n+i+1)+" and model "+objname)
        cmd.hide("everything","pp")
        cmd.delete("pp")
 

    def make_super_cell():
        global global_delocalized_bonds_list

        d = form.input_dimension.currentText()
        if d == "0":
           return 

        ratio_i = form.list_ratio.currentIndex()
        ratio = 0.1*(ratio_i+1)+0.01 

        make_unit_cell()

        #remove delocalized bonds in unit cell
        for l in range(len(global_delocalized_bonds_list)):
            i = global_delocalized_bonds_list[l][0]
            j = global_delocalized_bonds_list[l][1]
            cmd.delete("pka")
            cmd.delete("pkb")
            cmd.select("pka","id "+str(i+1)+" and "+"model geom")
            cmd.select("pkb","id "+str(j+1)+" and "+"model geom")
            cmd.unbond("pka","pkb")
            cmd.delete("pka")
            cmd.delete("pkb")

        model_o = cmd.get_model('geom', 1)
        elem = []
        for at in model_o.atom:
            elem.append(at.name)

        xyz_coor = cmd.get_model('geom', 1).get_coord_list()
        natom = len(xyz_coor)
        supercell_xyz = []
        supercell_xyz = append_xyz(supercell_xyz, xyz_coor)

        if int(d)>=1:    
           v1 = form.input_v1.text().split(",")
           v1 = [float(i) for i in v1]
           add_xyz = translate_geom( xyz_coor, v1, "+" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)

        if int(d)>=2:
           v2 = form.input_v2.text().split(",")
           v2 = [float(i) for i in v2]
           add_xyz = translate_geom( xyz_coor, v2, "+" ) 
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)

           add_xyz = translate_geom( translate_geom( xyz_coor, v1, "+" ), v2, "+" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)


        if int(d)>=3:
           v3 = form.input_v3.text().split(",")
           v3 = [float(i) for i in v3]
           add_xyz = translate_geom( xyz_coor, v3, "+" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)

           add_xyz = translate_geom( translate_geom( xyz_coor, v1, "+" ), v3, "+" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)

           add_xyz = translate_geom( translate_geom( xyz_coor, v2, "+" ), v3, "+" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)
           
           add_xyz = translate_geom( translate_geom( translate_geom( xyz_coor, v3, "+" ), v1, "+" ), v2, "+" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)



        write_xyz("supercell.xyz",supercell_xyz,elem)
        cmd.load("supercell.xyz","supercell")
        # delete file
        os.remove("supercell.xyz")

        set_style("supercell",1,"0.05","0.02") 

        cmd.enable("supercell")
        set_valence_obj("supercell",0)

        # show partial supercell
        MIN_o = get_min_corner(xyz_coor)
        #######cmd.pseudoatom(pos=MIN_o, object="basis_o")
        if int(d)>=1:
           for i in range(natom):
               p_xyz = xyz_coor[i]
               po = list_minus(p_xyz,MIN_o)
               ratio_x = list_dotp(po,v1)/list_dotp(v1,v1)
               if ratio_x > ratio:
                  hide_outer(1,natom,i,"supercell") 

        if int(d)>=2:
           for i in range(natom):
               p_xyz = xyz_coor[i]
               po = list_minus(p_xyz,MIN_o)
               ratio_x = list_dotp(po,v1)/list_dotp(v1,v1)
               ratio_y = list_dotp(po,v2)/list_dotp(v2,v2)
               if ratio_y > ratio:
                  hide_outer(2,natom,i,"supercell") 
               if ratio_x > ratio or ratio_y > ratio:
                  hide_outer(3,natom,i,"supercell") 

        if int(d)>=3:
           for i in range(natom):
               p_xyz = xyz_coor[i]
               po = list_minus(p_xyz,MIN_o)
               ratio_x = list_dotp(po,v1)/list_dotp(v1,v1)
               ratio_y = list_dotp(po,v2)/list_dotp(v2,v2)
               ratio_z = list_dotp(po,v3)/list_dotp(v3,v3)
               if ratio_z > ratio:
                  hide_outer(4,natom,i,"supercell") 
               if ratio_x > ratio or ratio_z > ratio:# 1,3
                  hide_outer(5,natom,i,"supercell") 
               if ratio_y > ratio or ratio_z > ratio:#2, 3
                  hide_outer(6,natom,i,"supercell") 
               if ratio_x > ratio or ratio_y > ratio or ratio_z > ratio:
                  hide_outer(7,natom,i,"supercell")

    def clicK_unit_cell():

        make_unit_cell()
        form.button_supercell.setEnabled(False)
        return 

    def make_unit_cell():
        global global_delocalized_bonds_list


        clear_cell()


        d = form.input_dimension.currentText()
        if d == "0":
           return         

        if d == "1":
          if "," in form.input_v1.text():
             v1 = form.input_v1.text().split(",")
          else:
             v1 = form.input_v1.text().split() 
          v1 = [float(i) for i in v1]
        if d == "2":
          if "," in form.input_v1.text():  
             v1 = form.input_v1.text().split(",")
             v2 = form.input_v2.text().split(",")
          else:
             v1 = form.input_v1.text().split() 
             v2 = form.input_v2.text().split()
          v1 = [float(i) for i in v1]
          v2 = [float(i) for i in v2]
        if d == "3":
          if "," in form.input_v1.text():
             v1 = form.input_v1.text().split(",")
             v2 = form.input_v2.text().split(",")
             v3 = form.input_v3.text().split(",")
          else:  
             v1 = form.input_v1.text().split()
             v2 = form.input_v2.text().split()
             v3 = form.input_v3.text().split()
          v1 = [float(i) for i in v1]
          v2 = [float(i) for i in v2]
          v3 = [float(i) for i in v3]

        xyz_coor = cmd.get_model('geom', 1).get_coord_list()
        MIN_Corner = get_min_corner(xyz_coor)

        cmd.pseudoatom(pos=MIN_Corner, object="basis_o")
        #print("created basis_o here...")
        #cmd.disable("basis_o")
        if int(d)>=1:
           cmd.pseudoatom(pos=list_add(MIN_Corner,v1), object="basis_1")
           cmd.distance("bv_1","basis_o","basis_1")
           basis_style("red","bv_1")
           #cmd.delete("basis_1") 

        if int(d)>=2:
           cmd.pseudoatom(pos=list_add(MIN_Corner,v2), object="basis_2")
           cmd.distance("bv_2","basis_o","basis_2")
           basis_style("green","bv_2")
           #cmd.delete("basis_2")

           cmd.pseudoatom(pos=list_add(list_add(MIN_Corner,v2),v1), object="basis_1p2")
           cmd.distance("bv_2_b","basis_1","basis_1p2")
           basis_style("green","bv_2_b")

           cmd.distance("bv_1_b","basis_2","basis_1p2")
           basis_style("red","bv_1_b")

        if int(d)==3:
           cmd.pseudoatom(pos=list_add(MIN_Corner,v3), object="basis_3")
           cmd.distance("bv_3","basis_o","basis_3")
           basis_style("blue","bv_3")

           cmd.pseudoatom(pos=list_add(v1,list_add(MIN_Corner,v3)), object="basis_3p1")
           cmd.distance("bv_3_b","basis_1","basis_3p1")
           basis_style("blue","bv_3_b")

           cmd.pseudoatom(pos=list_add(v2,list_add(MIN_Corner,v3)), object="basis_3p2")  
           cmd.distance("bv_3_c","basis_2","basis_3p2")
           basis_style("blue","bv_3_c")

           cmd.pseudoatom(pos=list_add(v1,list_add(v2,list_add(MIN_Corner,v3))), object="basis_3p1p2")
           cmd.distance("bv_3_d","basis_1p2","basis_3p1p2")
           basis_style("blue","bv_3_d")

           cmd.distance("bv_1_d","basis_3","basis_3p1")
           cmd.hide("label","bv_1_d")
           basis_style("red","bv_1_d")

           cmd.distance("bv_1_c","basis_3p2","basis_3p1p2")
           cmd.hide("label","bv_1_c")
           basis_style("red","bv_1_c")

           cmd.distance("bv_2_d","basis_3","basis_3p2")
           cmd.hide("label","bv_2_d")
           basis_style("green","bv_2_d")

           cmd.distance("bv_2_c","basis_3p1","basis_3p1p2")
           cmd.hide("label","bv_2_c")
           basis_style("green","bv_2_c")

        cmd.delete("basis_*") 

    def clear_cell():
        global global_delocalized_bonds_list
        cmd.delete("bv*")
        cmd.delete("supercell")
        cmd.delete("basis_*")
        cmd.enable("geom")      # added
        form.button_supercell.setEnabled(True)

    def start_anglewiz():
        anglewiz.reset()

        cmd.set_wizard(anglewiz)

    def remove_row():
        global global_row_index

        form.tableWidget.removeRow(global_row_index)
        form.button_remove.setEnabled(False)


    def select_row():
        global global_row_index 

        r = form.tableWidget.currentRow()
        global_row_index = r

        form.button_remove.setEnabled(True)
        return

    def title():

        logo="\
        o     o     o             8             .oo           o    o  \n\
        8     8b   d8             8            .P 8           8b   8  \n\
        8     8`b d'8 .oPYo. .oPYo8 .oPYo.    .P  8           8`b  8 .oPYo. odYo. .oPYo. \n\
        8     8 `o' 8 8    8 8    8 8oooo8   oPooo8   ooooo   8 `b 8 .oooo8 8' `8 8    8 \n\
        8     8     8 8    8 8    8 8.      .P    8           8  `b8 8    8 8   8 8    8 \n\
        8oooo 8     8 `YooP' `YooP' `Yooo' .P     8           8   `8 `YooP8 8   8 `YooP' \n\n\n"
                
        authorship = " "*39+" Ver 1.0.1 \n"+" "*38+" Feb 26, 2022 \n"+" "*28+" Developed by: Yunwen Tao, Ph.D.\n" 

        return logo,authorship

    def save_table():
        from datetime import datetime
        global global_save_text

        n = form.tableWidget.rowCount()
        if n == 0 :
           return
        maxl = 0
        for i in range(n):
           ll = len(form.tableWidget.item(i, 0).text())
           if maxl < ll:
              maxl = ll 
        filepath = browse_filename()
        #print(filepath)

        logo,author = title()
        
        f1 = open(filepath,"w")
        f1.write(logo)
        f1.write(author)



        # what to save
        # lmodea-nano title
        # 1. program type  [0]
        # 2. input file name [1], input file content [2]
        # 3. "geom" xyz [3]
        # 4. dimension [5]
        # 5. cell parameters [4] 
        # 6. table content
        # 7. eigenvalues [6] 
        f1.write("\n Program: "+global_save_text[0]+"\n")
        f1.write(" Input file: "+global_save_text[1]+"\n Input file content:\n")
        for i in range(len(global_save_text[2])):
            f1.write("  "+global_save_text[2][i]+"\n")

        f1.write("\n Dimension of the system: "+str(global_save_text[5])+"\n")
        f1.write("\n Lattice vectors (Angstrom):\n")
        for i in range(int(global_save_text[5])):
            f1.write("  "+global_save_text[4][i]+"\n")

        f1.write("\n Atomic positions (Angstrom):\n")   
        #print(global_save_text)
        for i in range(len(global_save_text[3])):
            f1.write("  "+global_save_text[3][i]+"\n")


        f1.write("\n First eigenvalues of Hessian:\n")
        for i in range(len(global_save_text[6])):
            f1.write("  "+str(global_save_text[6][i]))
        
        f1.write("\n\n")
       
        f1.write(" Local mode analysis result:\n")
        f1.write("  No. Atoms Length/Angle Force constant Frequency Comment    \n")
        for i in range(n): 
            f1.write('{0:5d}'.format(i+1))  
            q = form.tableWidget.item(i, 0).text()
            f1.write(" ")
            f1.write(q.ljust(maxl, ' '))
            qn= float(form.tableWidget.item(i, 1).text())
            f1.write('%13.3f' % qn)
            ka= float(form.tableWidget.item(i, 2).text())
            f1.write('%15.3f' % ka)
            w=int(form.tableWidget.item(i, 3).text())
            f1.write('{0:10d}'.format(w))
            
            if form.tableWidget.item(i, 4) is not None:
               cm = form.tableWidget.item(i, 4).text()
               f1.write(" "+cm)


            f1.write("\n")

        now = datetime.now()
        dt_string = now.strftime("%m/%d/%Y %H:%M:%S")
        f1.write("\n "+dt_string+"\n")

        f1.close()     

        return


    def start_bondwiz():
        bondwiz.reset()


        cmd.set_wizard(bondwiz)

        #

    def close_window():
        global global_delocalized_bonds_list
        global global_elem
        global global_fm
        global global_dim
        global global_mass
        global global_minv
        global global_save_text

        global_dim = 0
        global_delocalized_bonds_list = []
        global_elem = []
        global_fm = []
        global_mass = []
        global_minv = []
        global_save_text = []


        dialog.close()

    def linear_mol():
        global global_iflinmol
           
        if form.check_linear.isChecked():
           global_iflinmol = True
           #print("linear mol")
        else:   
           global_iflinmol = False
           #print("non-linear mol")



    def about_window():
        ret = QtWidgets.QMessageBox.question(None,'About LModeA-nano', "Author: Yunwen Tao, Ph.D.\nEmail: ywtao.smu@gmail.com", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,)
        if ret == QtWidgets.QMessageBox.Yes:
           print("Have a good day!")


    form.button_save.clicked.connect(save_table)
    form.button_remove.clicked.connect(remove_row) 
    form.tableWidget.itemClicked.connect(select_row)
    form.button_about.clicked.connect(about_window)
    form.button_angle.clicked.connect(start_anglewiz)
    form.button_bond.clicked.connect(start_bondwiz)
    form.button_supercell.clicked.connect(make_super_cell)
    form.button_clearCell.clicked.connect(clear_cell)
    form.button_unitCell.clicked.connect(clicK_unit_cell)#(make_unit_cell)
    form.radio_mol.toggled.connect(select_mol) 
    form.button_confirm.clicked.connect(update_dimension) 
    form.button_load.clicked.connect(loadxyz)#_temp)
    form.button_find.clicked.connect(find_input_file)
    form.button_close.clicked.connect(close_window)#dialog.close)
    form.check_linear.stateChanged.connect(linear_mol)

    return dialog


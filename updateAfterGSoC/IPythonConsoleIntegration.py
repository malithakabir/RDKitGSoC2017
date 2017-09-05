
"""RDKit Conformer Browser
Derived by Malitha Humayun Kabir from Greg Landrum's code as a part of GSoC 2017
Project : RDKit - 3Dmol.js integration
Mentors: Paul Czodrowski and Greg Landrum
Acknowledgement: Peter Gedeck reviewed, provided advice on restructuring, and wrote initial MolViewState class.
Date: 5th September 2017
Email# malitha12345@gmail.com
"""

try:
    import py3Dmol
    _canUse3D = True
except ImportError:
    _canUse3D = False
    
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from ipywidgets import Layout, Label, Button, Box, HBox, VBox
from ipywidgets import Dropdown, HTML, Checkbox, Button
from IPython.display import display
from IPython.display import HTML as scriptHTML
import time
import sys
import copy

PY3 = sys.version_info[0] >= 3

BGCOLORS_3D = ('0x000000', '0xeeeeee', '0xffffff')

PROP_RDKIT = tuple(sorted(prop for prop, _ in Descriptors._descList))

DRAWING_LIGAND_3D=('line', 'cross', 'stick', 'sphere', 'surface', 'ballstick')

DRAWING_PROTEIN_3D=('line', 'cartoon', 'surface', 'cartoonWithTube')

COLOR_SCHEME_3D=('default', 'greenCarbon', 'cyanCarbon', 'magentaCarbon',
                 'yellowCarbon', 'whiteCarbon', 'orangeCarbon', 'purpleCarbon', 
                 'blueCarbon', 'ssPyMOL', 'ssJmol', 'Jmol', 'amino', 
                 'shapely', 'nucleic', 'chain', 'chainHetatm', 'prop')
# TODO: implement COLOR_3D
COLOR_3D = ('white', 'silver', 'gray', 'black', 'red', 'maroon', 'yellow', 
            'orange', 'olive', 'lime', 'green', 'aqua', 'cyan', 'teal', 'blue', 
            'navy', 'fuchsia', 'magenta', 'purple')

def Check3Dmolpromise():
    """This function checks whether 3Dmol.js loaded or not."""
    uniqueid = str(time.time()).replace('.','')
    input_form = """<div id="%s"></div>""" % uniqueid
    javascript = """
    <script type="text/Javascript">
        function checkLibraryLoad()
            {
                if(typeof $3Dmolpromise === 'undefined')
                {
                    res = '3DMol.js not loaded'
                } else {
                    res = '3DMol.js already loaded'
                }
                return res
            }
        document.getElementById('%s').innerHTML = checkLibraryLoad()
    </script>
    """ % uniqueid
    return(scriptHTML(input_form + javascript))
    
def ProcessLigandDict(ligandDict, createChildKey = 'parent'):
    """This function adds another key to the dictionary of the molecule."""
    if isinstance(ligandDict, dict) is False:
        raise TypeError("ligandDict must be a dictionary")
        
    keys = list(ligandDict.keys())
    firstKey = keys[0]
    
    if isinstance(ligandDict[firstKey], dict):
        raise TypeError("ProcessLigandDict doesn't support nested ligandDict")
        
    newLigandDict = {}
    
    for molId in keys:
        newLigandDict[molId] = {}
        newLigandDict[molId][createChildKey] = ligandDict[molId]
    return newLigandDict
    
def CopyMolWithChildKey(ligandDict, copyFrom = 'parent', newChildKey = 'minimized'):
    """This function adds another key to the dictionary of the molecule."""
    if isinstance(ligandDict, dict) is False:
        raise TypeError("ligandDict must be a dictionary")
        
    keys = list(ligandDict.keys())
    
    for molId in keys:
        oldMol = ligandDict[molId][copyFrom]
        newMol = copy.deepcopy(oldMol)
        ligandDict[molId][newChildKey] = newMol
        
def ConfGenInLigandDict(ligandDict, onChildKey = 'minimized', numConfs = 10, molIds = 'allMols'):
    """This function adds another key to the dictionary of the molecule."""
    if molIds is 'allMols':
        molIds = list(ligandDict.keys())
        
    params = AllChem.ETKDG()
    params.numThreads=3
    for molId in molIds:
        AllChem.EmbedMultipleConfs(ligandDict[molId][onChildKey], numConfs=numConfs, params=params)
        
def MinimizeLigand(ligandDict, onChildKey = 'minimized',
                   createPerConfDataKey = 'energy',
                   molAndConfIds = 'allConfs', 
                   ff = 'UFF', 
                   maxIters = 50):
    """ This function takes a dictionary of ligand and does energy minimization to the ligands"""
    if ff not in ('MMFF', 'UFF'):
        raise TypeError("ff can be either MMFF or UFF")
        
    if isinstance(ligandDict, dict) is False:
        raise TypeError("ligandDict must be dict")
        
    if all(isinstance(key, str if PY3 else basestring) for key in ligandDict.keys()) is False:
        raise TypeError("keys of ligandDict must be str")
        
    if molAndConfIds is 'allConfs':
        molAndConfIdsTuple = ()
        for molId in list(ligandDict.keys()):
            mol = ligandDict[molId][onChildKey]
            confIds = list(range(mol.GetNumConformers()))
            for cid in confIds:
                molAndConfIdsTuple = molAndConfIdsTuple + ((molId, cid),)
    else:
        molAndConfIdsTuple = molAndConfIds
        
    if all(isinstance(molId, str if PY3 else basestring) for (molId, confId) in molAndConfIdsTuple) is False:
        raise TypeError("keys of ligandDict must be str")
        
    for Ids in molAndConfIdsTuple:
        
        molId = Ids[0]
        confId = Ids[1]
        
        energy = ligandDict[molId].setdefault(createPerConfDataKey, {})
        mol = ligandDict[molId][onChildKey]
        
        if ff == 'MMFF':
            getFF = AllChem.MMFFGetMoleculeForceField(mol, confId=confId)
            getFF.Minimize(maxIts = maxIters)
            energy[confId] = getFF.CalcEnergy()
            
        elif ff == 'UFF':
            getFF = AllChem.UFFGetMoleculeForceField(mol, confId=confId)
            getFF.Minimize(maxIts = maxIters)
            energy[confId] = getFF.CalcEnergy()
        else:
            raise ValueError("invalid conformer id")
            
def AddPropToLigandDict(ligandDict, onChildKey = 'minimized'):
    """ Add property to the mol """
    for molId in ligandDict:
        mol = ligandDict[molId][onChildKey]
        for prop_name in PROP_RDKIT:
            calculator = Descriptors.__dict__[prop_name]
            mol.SetProp(prop_name, str(calculator(mol)))
            
class MolViewState(object):
    def __init__(self, 
                 molecules, 
                 protein, 
                 additionalChildKeysForMolRender = ['parent'],
                 childKeyForConfSelection = 'minimized',
                 childKeyForDataPerConf = 'energy',
                 childKeyForPropSelection = 'minimized'):
        """ Molecules is dictionary of molecules """
        
        self.ligandDict = molecules
        self.protein = protein
        
        self.additionalChildKeysForMolRender = additionalChildKeysForMolRender
        self.childKeyForConfSelection = childKeyForConfSelection
        
        self.childKeyForDataPerConf = childKeyForDataPerConf
        self.childKeyForPropSelection = childKeyForPropSelection
        # These should have reasonable initial values
        self.rdkit_mol_select = set()
        self.rdkit_conf_select = set()
        
    def selectAllMols(self):
        """ Select all molecules"""
        self.rdkit_mol_select = set(self.ligandDict)
        
    def selectAllConfs(self):
        """ For all selected molecules, select all conformations"""
        for molIds in self.selectedMolNames:
            
            mol = self.ligandDict[molIds][self.childKeyForConfSelection]
            nConformers = mol.GetNumConformers()
            
            if nConformers > 1:
                self.rdkit_conf_select = set(range(nConformers))
            else:
                self.rdkit_conf_select = {0}
                
    def selectSingleMol(self, molId):
        """ Select molecule"""
        self.rdkit_mol_select = {molId}
        
    def selectSingleConf(self, confId):
        """ Select conformation"""
        self.rdkit_conf_select = {confId}
        
    def generateIds(self):
        """generate paired id tuple from rdkit_mol_select and rdkit_conf_select"""
        self.idPaired = set()
        for molId in self.rdkit_mol_select:
            for confId in self.rdkit_conf_select:
                self.idPaired.add((molId,confId))
                
    def update_rdkit_select_unique_ids(self):
        """generate unique molId and confId from paired id tuple"""
        self.rdkit_mol_select = set()
        self.rdkit_conf_select = set()
        for (molId, confId) in self.idPaired:
            self.rdkit_mol_select.add(molId)
            self.rdkit_conf_select.add(confId)
            
    def deleteIds(self, molAndConfIds):
        """ takes directly parsed paired id and store """
        for ids in molAndConfIds:
            self.idPaired.discard(ids)
        self.update_rdkit_select_unique_ids()
        
    def addIds(self, molAndConfIds):
        """ takes directly parsed paired id and store """
        for ids in molAndConfIds:
            self.idPaired.add(ids)
        self.update_rdkit_select_unique_ids()
        
    def selectSingleId(self, molAndConfIds):
        """ takes directly parsed paired id and store """
        self.idPaired = set()
        for ids in molAndConfIds:
            self.idPaired.add(ids)
        self.update_rdkit_select_unique_ids()
        
    @property
    def selectedModels(self):
        """ Iterator over all selected models (molecules/conformations) """
        modelDict = {}
        k_for_conf = self.childKeyForConfSelection
        k_for_additional = self.additionalChildKeysForMolRender
        
        for (molId, confId) in self.idPaired:
            molData = self.ligandDict[molId]
            confKey = k_for_conf in list(molData.keys())
            if confKey:
                mol_with_conf = molData[k_for_conf]
                modelDict[k_for_conf] = Chem.MolToMolBlock(mol_with_conf,confId=confId)
                for key in k_for_additional:
                    additionalkey = key in list(molData.keys())
                    if additionalkey:
                        mol_additional = molData[key]
                        if confId in list(range(mol_additional.GetNumConformers())):
                            modelDict[key] = Chem.MolToMolBlock(mol_additional,confId=confId)
            yield modelDict
            
    @property
    def selectedMolNames(self):
        """ Return the names of all selected molecules """
        return self.rdkit_mol_select
        
    @property
    def selectedConfIds(self):
        """ Return the names of all selected confIds """
        return self.rdkit_conf_select
        
    @property
    def selectedMolecules(self):
        """ Return the selected molecules """
        return [self.ligandDict[name][self.childKeyForConfSelection] for name in self.selectedMolNames]
        
    @property
    def allConfIds(self):
        """ Return the number of conformations - use the first selected molecule to determine """
        nconfIds = self.selectedMolecules[0].GetNumConformers()
        return list(range(nconfIds))
        
    @property
    def getPropPrecalculated(self):
        """ Return the precalculated properties """
        if len(self.selectedMolNames) != 1:
            return None
        else:
            mol = self.ligandDict[list(self.selectedMolNames)[0]][self.childKeyForPropSelection]
            
            if len(mol.GetPropNames()) == 0:
                return 'Not found'
            else:
                return {prop: mol.GetProp(prop) for prop in mol.GetPropNames()}
                
    @property
    def getMinimizationEnergy(self):
        """ Return the minimization data """
        if len(self.selectedMolNames) != 1:
            return None
        else:
            molData = self.ligandDict[list(self.selectedMolNames)[0]]
            if self.childKeyForDataPerConf in molData:
                return molData[self.childKeyForDataPerConf]
            else:
                return 'childKeyForDataPerConf not found'
                
    @property
    def atomLabel(self):
        """ Return the atomLabel """
        if len(self.selectedMolNames) == 1 and len(self.selectedConfIds) == 1:
            
            mol = self.ligandDict[list(self.selectedMolNames)[0]][self.childKeyForConfSelection]
            sconf = mol.GetConformer(list(self.selectedConfIds)[0])
            xyz = sconf.GetPositions()
            atomLabel = []
            for i in list(range(sconf.GetNumAtoms())):
                label = ('{}{}'.format(mol.GetAtomWithIdx(i).GetSymbol(), i), 
                         xyz[i]
                        )
                atomLabel.append(label)
            return atomLabel
    
    @property
    def confLabel(self):
        """ Return the conformer label """
        if len(self.selectedMolNames) == 1 and len(self.selectedConfIds) == 1:
            confLabel = list(self.selectedMolNames)[0] + ':' + str(list(self.selectedConfIds)[0])
            return confLabel
            
            
        
        
class MolView3D(object):
    
    def __init__(self, 
                 ligandDict = None,
                 protein = None, 
                 childKeysForMolRender = ['parent', 'minimized'], 
                 childKeyForConfSelection = 'minimized',
                 childKeyForDataPerConf = 'energy',
                 childKeyForPropSelection = 'minimized',
                 ligStyle = ['stick', 'stick'], protStyle='cartoon',
                 ligColor = ['default','orangeCarbon'], protColor='ssJmol',
                 molVisibilityPanel = True, 
                 stylePanel = None, 
                 labelPanel = False,
                 propertyPanel = False,
                 emPanel = False,
                 molAndConfIds = 'allConfs'):
        """This function initiates required widgets and 3Dmol.js viewer"""
        if stylePanel is not None:
            if stylePanel not in ('lig', 'prot', 'ligprot'):
                raise KeyError('stylePanel can be: lig, prot, ligprot')
                
        if type(childKeysForMolRender) is not list:
            raise TypeError('childKeysForMolRender must be a list')
            
        self.wgNameForLigandVisible = []
        self.wgNameForLigandStyle = []
        self.wgNameForLigandColor = []
        self.onStart = True
        self.wg_values = {}
        
        
        self.childKeysForMolRender = childKeysForMolRender
        self.childKeyForConfSelection = childKeyForConfSelection
        
        self.childKeyForDataPerConf = childKeyForDataPerConf
        self.childKeyForPropSelection = childKeyForPropSelection
        
        self.ligStyleScheme = {}
        self.ligColorScheme = {}
        for (ids,lig) in enumerate(childKeysForMolRender):
            self.ligStyleScheme[lig] = ligStyle[ids]
            self.ligColorScheme[lig] = ligColor[ids]
            
        self.protStyle = protStyle
        self.protColorScheme = protColor
        
        self.molVisibilityPanel = molVisibilityPanel
        self.propertyPanel = propertyPanel
        self.stylePanel = stylePanel
        self.labelPanel = labelPanel
        self.emPanel = emPanel
        
        
        # Right hand panel (widgets)
        self.style = {'description_width': 'initial'}
        self.rightBoxWidth = '45%'
        self.wg_height = '30px'
        self.wgBox = []
        self.itemLayout=Layout(display='flex', flex_flow='row', 
                               justify_content='flex-start', 
                               height=self.wg_height)
        self.widgetsFor3DView = []
        
        self.allWidgets = ['wg_selectedMolsConfsView',
                           'wg_selectConfIdfrom', 
                           'wg_selectAllMols', 'wg_selectAllConfs',
                           'wg_molId', 'wg_confId', 
                           'wg_proteinVisible',
                           'wg_propSelect', 'wg_propView', 
                           'wg_energyView', 
                           'wg_ligStyle', 'wg_ligColorScheme',
                           'wg_protStyle', 'wg_protColorScheme',
                           'wg_confLabel', 'wg_atomLabel', 
                           'wg_background']
        
        # Molecule selection panel
        self.MoleculeSelectionPanelUI()
        
        # Molecule visibility (show/hide)
        if self.molVisibilityPanel:
            self.MoleculeVisibilityPanelUI(ligandDict, protein)
            
        # property show
        if propertyPanel:
            self.PropertyViewPanelUI()
        
        # minimization data show
        if self.emPanel:
            self.EnergyMinimizationDataViewPanelUI()
            
        # stylePanel
        if stylePanel is not None:
            if stylePanel:
                self.StylePanelUI()
        
        # add conf labels and atom labels 
        if labelPanel:
            self.LabelPanelUI()
        
        
        # background color
        self.BackgroundColorPanelUI()
        
        # Add observer to interactive widgets
        self.AddObserverToWidgets()
        
        # Rendering left (molecule viewer) and right (all widgets) panels
        self.RenderWidgetsWithViewer()
        
        
        # adding model to viewer
        self.SetMolData(ligandDict, protein, 
                        childKeysForMolRender, 
                        childKeyForConfSelection, 
                        childKeyForDataPerConf,
                        childKeyForPropSelection,
                        molAndConfIds = None)
        
    def GetWidgetValue(self):
        """this captures values from all the widgets except wg_modelSelect, wg_modelAdd and wg_zoomTo"""
        for wg in self.allWidgets:
            if wg in self.__dict__:
                self.wg_values[wg] = self.__dict__[wg].value
                
    def MoleculeSelectionPanelUI(self):
        """ligand selection widgets"""
        self.wg_selectedMolsConfsView = HTML(description='', value='')
        b1 = Box([Label(value='ids'),self.wg_selectedMolsConfsView], layout=self.itemLayout)
        
        tempLayout = Layout(display='flex', justify_content='flex-start', width='150px', height=self.wg_height)
        
        self.wg_selectConfIdfrom = Dropdown(description='', options=self.childKeysForMolRender,
                                            value=self.childKeyForConfSelection, layout=tempLayout)
        self.wg_updateConfId = Button(description="updateConfId", button_style='success',
                                      layout=Layout(width='100px'))
        box_for_confId_selection = HBox([self.wg_selectConfIdfrom, self.wg_updateConfId])
        b2 = Box([Label(value='confId from'),box_for_confId_selection], layout=self.itemLayout)
        
        self.wg_selectAllMols = Checkbox(description='selectAllMols', value=False)
        self.wg_selectAllConfs = Checkbox(description='selectAllConfs', value=False)
        box_for_allmol_allconf = HBox([self.wg_selectAllMols, self.wg_selectAllConfs])
        b3 = Box([Label(value=''), box_for_allmol_allconf], layout=self.itemLayout)
        
        molIdLayout = Layout(display='flex', justify_content='flex-start', width='200px', height=self.wg_height)
        confIdLayout = Layout(display='flex', justify_content='flex-start', width='140px', height=self.wg_height)
        
        self.wg_molId = Dropdown(description='molId', options=['select'],value='select', layout=molIdLayout)
        self.wg_confId = Dropdown(description='confId', options=['select'], value='select', layout=confIdLayout)
        b4 = Box([HBox([self.wg_molId, self.wg_confId])],layout=self.itemLayout)
        
        self.wg_modelDelete = Button(description='delete', button_style='success', layout=Layout(width='100px'))
        self.wg_modelAdd = Button(description='add', button_style='success', layout=Layout(width='100px'))
        self.wg_modelSelect = Button(description='select', button_style='success', layout=Layout(width='100px'))
        b5 = Box([Label(value=''),HBox([self.wg_modelDelete, self.wg_modelAdd, self.wg_modelSelect])
                 ],layout=self.itemLayout)
        
        self.wg_zoomTo = Button(description="zoomTo", button_style='success', layout=Layout(width='150px'))
        self.wg_reDraw = Button(description='reDraw', button_style='success', layout=Layout(width='150px'))
        b6 = Box([Label(value=''),HBox([self.wg_zoomTo, self.wg_reDraw])],layout=self.itemLayout)
        
        molConfBox = VBox([b1, b2, b3, b4, b5, b6], layout=Layout(border='solid', border_color = 'blue'))
        
        self.wgBox.append(molConfBox)
        
    def MoleculeVisibilityPanelUI(self, ligandDict, protein):
        """ligand, protein, and energy minimized ligand visibility widgets"""
        protein_sh = True if protein is not None else False
        if protein_sh:
            self.wg_proteinVisible = Checkbox(description='proteinVisible', value=True)
            self.wgBox.append(Box([Label(value=''),self.wg_proteinVisible], layout=self.itemLayout))
            
        ligand_sh = True if ligandDict is not None else False    
        if ligand_sh:
            for name in self.childKeysForMolRender:
                wgName = 'wg_' + name + 'Visible'
                self.__dict__[wgName] = Checkbox(description = wgName[3:], value=True)
                self.wgBox.append(Box([Label(value=''), self.__dict__[wgName]], layout=self.itemLayout))
                self.wgNameForLigandVisible.append(wgName)
                self.allWidgets.append(wgName)
                    
    def PropertyViewPanelUI(self):
        """widgets for ligand property selection and property show"""
        self.wg_propView = HTML(description='prop', value='initializing...', style = self.style)
        self.wgBox.append(Box([Label(value=''),self.wg_propView], layout=self.itemLayout))
        
        self.wg_propSelect = Dropdown(description='prop',options=['select'], value='select', style = self.style)
        self.wgBox.append(Box([Label(value=''),self.wg_propSelect], layout=self.itemLayout))
        
    def EnergyMinimizationDataViewPanelUI(self):
        """widgets for showing energy data for each conformer"""
        self.wg_energyView = HTML(description='', value='initializing...')
        self.wgBox.append(Box([Label(value=''),self.wg_energyView], layout=self.itemLayout))
        
    def LigStylePanelUI(self):
        """widgets for selecting ligand representation and coloring scheme"""
        for (ids,name) in enumerate(self.childKeysForMolRender):
            
            wgNameForStyle = 'wg_' + name + 'Style'
            self.wgNameForLigandStyle.append(wgNameForStyle)
            
            self.__dict__[wgNameForStyle] = Dropdown(description=wgNameForStyle[3:], 
                                                     options=DRAWING_LIGAND_3D,
                                                     value=self.ligStyleScheme[name], style = self.style)
            self.wgBox.append(Box([Label(value=''), 
                                   self.__dict__[wgNameForStyle]
                                  ], layout=self.itemLayout))
            
            wgNameForColor = 'wg_' + name + 'Color'
            self.wgNameForLigandColor.append(wgNameForColor)
            
            self.__dict__[wgNameForColor] = Dropdown(description=wgNameForColor[3:], 
                                                     options=COLOR_SCHEME_3D,
                                                     value=self.ligColorScheme[name], style = self.style)
            self.wgBox.append(Box([Label(value=''), 
                                   self.__dict__[wgNameForColor]
                                  ], layout=self.itemLayout))
            self.allWidgets.extend([wgNameForStyle, wgNameForColor])
            
    def ProteinStylePanelUI(self):
        """widgets for selecting protein representation"""
        self.wg_protStyle = Dropdown(description='protStyle', 
                                     options=DRAWING_PROTEIN_3D, 
                                     value=self.protStyle,
                                     style = self.style)
        self.wgBox.append(Box([Label(value=''),self.wg_protStyle], layout=self.itemLayout))
        
        self.wg_protColorScheme = Dropdown(description='protColor', 
                                           options=COLOR_SCHEME_3D, 
                                           value=self.protColorScheme, 
                                           style = self.style)
        self.wgBox.append(Box([Label(value=''),self.wg_protColorScheme], layout=self.itemLayout))
        
    def StylePanelUI(self):
        """widgets for selecting which style panel will be appeared"""
        if self.stylePanel == 'lig':
            self.LigStylePanelUI()
        if self.stylePanel == 'prot':
            self.ProteinStylePanelUI()
        if self.stylePanel == 'ligprot':
            self.LigStylePanelUI()
            self.ProteinStylePanelUI()
            
    def LabelPanelUI(self):
        """widgets for labeling ligand (not energy minimized ligand)"""
        self.wg_confLabel = Checkbox(description='confLabel', value=False)
        self.wg_atomLabel = Checkbox(description='atomLabel', value=False)
        box_for_labeling = HBox([self.wg_confLabel,self.wg_atomLabel])
        self.wgBox.append(Box([Label(value=''), box_for_labeling], layout=self.itemLayout))
        
    def BackgroundColorPanelUI(self):
        """widgets for selecting viewer background"""
        self.wg_background = Dropdown(description='background', 
                                      options=BGCOLORS_3D, 
                                      value=BGCOLORS_3D[0], style = self.style)
        self.wgBox.append(Box([Label(value=''),self.wg_background], layout=self.itemLayout))
        
    def AddObserverToWidgets(self):
        """This function sets observer to all the interactive widgets except zoomTo"""
        self.wg_updateConfId.on_click(self.handle_updateConfId_button)
        self.wg_zoomTo.on_click(self.handle_zoomTo_button)
        self.wg_reDraw.on_click(self.handle_reDraw_button)
        self.wg_modelDelete.on_click(self.handle_modelDelete_button)
        self.wg_modelAdd.on_click(self.handle_modelAdd_button)
        self.wg_modelSelect.on_click(self.handle_modelSelect_button)
        self.wg_selectAllMols.observe(self.handle_selectAllMols_checkbox, names='value')
        self.wg_selectAllConfs.observe(self.handle_selectAllConfs_checkbox, names='value')
        
    def RenderWidgetsWithViewer(self):
        """this creates a html table having an id as and then insert the 3DMol.js viewer in table"""
        tableUID = str(time.time()).replace('.','')
        
        # left panel (container holding table for viewer)
        size = (435, 485)
        viewerLayout=Layout(width=str(size[0]+5)+'px', height=str(size[1]+5)+'px', border='solid')
        wg_leftWidget = HTML('''<table><tr><td id="%s"></td></tr></table>'''% tableUID, layout=viewerLayout)
        
        # right panel
        wg_rightBox=VBox(self.wgBox, layout=Layout(border='solid', width=self.rightBoxWidth))
        
        # Combining left and right panels
        viewer = HBox([wg_leftWidget, wg_rightBox])
        
        # displaying everything
        display(viewer)
        
        # inserting 3DMol.js viewer in existing container (table)
        self.view = py3Dmol.view(width=size[0],height=size[1])
        self.view.setBackgroundColor('0x000000')
        self.view.zoomTo()
        
        display(self.view.insert(tableUID))
        
    def handle_updateConfId_button(self, b):
        """This function handles modelDelete button"""
        self.UpdateConfId(childKeyForConfSelection = self.wg_selectConfIdfrom.value)
        
    def handle_reDraw_button(self, b):
        """This function handles reDraw button"""
        self.render3D()
        
    def handle_zoomTo_button(self, b):
        """This function handles zoomTo button"""
        self.view.zoomTo()
        display(self.view.update())
        
    def handle_modelDelete_button(self, b):
        """This function handles modelDelete button"""
        molAndConfIds = ((self.wg_molId.value, self.wg_confId.value),)
        self.DeleteMolAndConf(molAndConfIds = molAndConfIds)
        
    def handle_modelAdd_button(self, b):
        """This function handles modelAdd button"""
        molAndConfIds = ((self.wg_molId.value, self.wg_confId.value),)
        self.AddMolAndConf(molAndConfIds = molAndConfIds)
        
    def handle_modelSelect_button(self, b):
        """This function handles modelSelect button"""
        molAndConfIds = ((self.wg_molId.value, self.wg_confId.value),)
        self.SelectMolAndConf(molAndConfIds = molAndConfIds)
    
    def handle_selectAllMols_checkbox(self, change):
        """This function handles selectAllMols checkbox"""
        self.SelectAllMols(confId = self.wg_confId.value, callFrom = 'widget')
        
    def handle_selectAllConfs_checkbox(self, change):
        """This function handles selectAllConfs checkbox"""
        self.SelectAllConfs(molId = self.wg_molId.value, callFrom = 'widget')
        
    def reInstantiateViewer(self):
        """This function resets mol data"""
        self.SetMolData(ligandDict = self.molViewState.ligandDict,
                        protein = self.molViewState.protein, 
                        childKeysForMolRender  = self.childKeysForMolRender ,
                        childKeyForConfSelection = self.childKeyForConfSelection,
                        childKeyForDataPerConf = self.childKeyForDataPerConf,
                        childKeyForPropSelection = self.childKeyForPropSelection,
                        molAndConfIds = None, reinstantiated = True)
                        
    def UpdateConfId(self, childKeyForConfSelection = 'minimized'):
        """This function update the childKeys from which conf ids rendered"""
        self.childKeyForConfSelection = childKeyForConfSelection
        self.wg_selectConfIdfrom.value = childKeyForConfSelection
        self.reInstantiateViewer()
        
    def updateWidgets(self, callFrom):
        """ update mol and conf id viewer"""
        if callFrom == 'SelectAllMols':
            self.molViewState.generateIds()
            
        if callFrom == 'SelectAllConfs':
            self.molViewState.generateIds()
            
        #For multiple molecules with unequal number of conformers, it is hard to determine acceptable confIds
        self.wg_confId.options = self.molViewState.allConfIds
        self.wg_selectedMolsConfsView.value = ', '.join([str(x) for x in self.molViewState.idPaired])
        if callFrom == 'SelectMolAndConf':
            if self.onStart:
                self.render3D()
                
    def SelectAllMols(self, confId = None, callFrom = 'nonWidget'):
        """ instantiates mol and conformer selection function of the MolViewState"""
        self.molViewState.selectAllMols()
        if confId is not None:
            self.molViewState.selectSingleConf(confId)
        if callFrom == 'widget':
            self.updateWidgets(callFrom = 'SelectAllMols')
            
    def SelectAllConfs(self, molId = None, callFrom = 'nonWidget'):
        """ instantiates mol and conformer selection function of the MolViewState"""
        if molId is not None:
            self.molViewState.selectSingleMol(molId)
        self.molViewState.selectAllConfs()
        if callFrom == 'widget':
            self.updateWidgets(callFrom = 'SelectAllConfs')
            
    def DeleteMolAndConf(self, molAndConfIds):
        """ instantiates mol and conformer selection function of the MolViewState"""
        self.molViewState.deleteIds(molAndConfIds = molAndConfIds)
        self.updateWidgets(callFrom = 'DeleteMolAndConf')
        
    def AddMolAndConf(self, molAndConfIds):
        """ instantiates mol and conformer selection function of the MolViewState"""
        self.molViewState.addIds(molAndConfIds = molAndConfIds)
        self.updateWidgets(callFrom = 'AddMolAndConf')
        
    def SelectMolAndConf(self, molAndConfIds):
        """ instantiates mol and conformer selection function of the MolViewState"""
        self.molViewState.selectSingleId(molAndConfIds = molAndConfIds)
        self.updateWidgets(callFrom = 'SelectMolAndConf')
        
    def SetMolData(self, 
                   ligandDict = None, 
                   protein = None, 
                   childKeysForMolRender = ['parent', 'minimized'], 
                   childKeyForConfSelection = 'minimized',
                   childKeyForDataPerConf = 'energy',
                   childKeyForPropSelection = 'minimized',
                   molAndConfIds = None, reinstantiated = False):
        """This function sets ligand dictionary, protein, and dict keys and initiates MolViewState class"""
        
        if isinstance(ligandDict, dict) is False:
            raise TypeError("ligandDict must be dict")
            
        if all(isinstance(key, str if PY3 else basestring) for key in ligandDict.keys()) is False:
            raise TypeError("keys of ligandDict must be str")
            
        self.childKeysForMolRender = childKeysForMolRender
        self.childKeyForConfSelection = childKeyForConfSelection
        
        additionalChildKeysForMolRender = set(self.childKeysForMolRender)
        additionalChildKeysForMolRender.discard(self.childKeyForConfSelection)
        self.additionalChildKeysForMolRender = list(additionalChildKeysForMolRender)
        
        self.childKeyForDataPerConf = childKeyForDataPerConf
        self.childKeyForPropSelection = childKeyForPropSelection
        
        if reinstantiated:
            self.molViewState = None
            
        self.molViewState = MolViewState(ligandDict, protein, 
                                         self.additionalChildKeysForMolRender, 
                                         childKeyForConfSelection, 
                                         childKeyForDataPerConf,
                                         childKeyForPropSelection)
        if self.molViewState.ligandDict is not None:
            
            keys = list(self.molViewState.ligandDict.keys())
            self.onStartIds = ((keys[0], 0),)
            
            self.wg_molId.options = keys
            self.wg_molId.value = keys[0]
            self.SelectMolAndConf(molAndConfIds = self.onStartIds)
            
    def ShowLigandProperty(self):
        """ Handles property in render3D function """
        preCalcProp = self.molViewState.getPropPrecalculated
        if isinstance(preCalcProp, dict):
            self.wg_propSelect.options = sorted(list(preCalcProp.keys()))
            prop=self.wg_propSelect.value
            self.wg_propView.value = prop + ' : ' + str(preCalcProp[prop])
        elif preCalcProp == 'Not found':
            self.wg_propView.value = 'data not found!'
        else:
            self.wg_propView.value = 'single molecule selection required!'
            
    def ShowMinimizationEnergy(self):
        """ Handles energy minimization data in render3D function """
        energy = self.molViewState.getMinimizationEnergy
        energyDataKey = self.molViewState.childKeyForDataPerConf
        confId = self.wg_confId.value
        
        if isinstance(energy, dict) and confId in energy:
            self.wg_energyView.value = 'em data ('+ energyDataKey +') : ' + str(energy[confId])
        elif isinstance(energy, dict) and confId not in energy:
            self.wg_energyView.value = 'data not found'
        elif energy is None:
            self.wg_energyView.value = 'single molecule selection required'
        else:
            self.wg_energyView.value = 'data not found'
            
    def AddLigandLabels(self):
        """ Handles ligand labeling (conf label and atom label) in render3D function """
        
        if self.wg_confLabel.value:
            if self.molViewState.confLabel is not None:
                self.view.addLabel(self.molViewState.confLabel, {'backgroundColor':'gray', 
                                                                 'fontColor':'white',
                                                                 'showBackground':'true',
                                                                 'alignment':'bottomCenter'})
        if self.wg_atomLabel.value:
            if self.molViewState.atomLabel is not None:
                for (label,xyz) in self.molViewState.atomLabel:
                    self.view.addLabel(label, {'inFront' : 'false', 
                                               'fontSize' : '12',
                                               'fontColor':'gray',
                                               'showBackground':'false',
                                               'position' : {'x':xyz[0], 'y':xyz[1], 'z':xyz[2]}
                                              })
                    
    def AddLigandStyle(self, name, modelId):
        """ Handles ligand and energy minimized ligand drawing style and color in AddLigandWithStyle function """
        wgNameStyle = 'wg_' + name + 'Style'
        wgNameColor = 'wg_' + name + 'Color'
        
        style = self.__dict__[wgNameStyle].value if wgNameStyle in self.__dict__  else self.ligStyleScheme[name]
        color = self.__dict__[wgNameColor].value if wgNameColor in self.__dict__ else self.ligColorScheme[name]
        
        coloring = 'colorscheme'
        if style == 'surface':
            self.view.addSurface('SES', {'model':modelId});
        elif style == 'ballstick':
            self.view.setStyle({'model':modelId},{'stick':{'radius':'0.2', coloring: color},
                                                  'sphere':{'radius':'0.4', coloring: color}
                                                 });
        else:
            self.view.setStyle({'model': modelId},{style:{coloring: color}})
            
    def AddProteinStyle(self):
        """ Handles protein drawing style in AddProteinWithStyle function """
        style = self.wg_protStyle.value if 'wg_protStyle' in self.__dict__ else self.protStyle
        color = self.wg_protColorScheme.value if 'wg_protColorScheme' in self.__dict__ else self.protColorScheme
        coloring = 'colorscheme'
        
        if style == 'surface':
            self.view.addSurface('SES', {'model': self.proteinModelId, coloring: color});
            
        elif style == 'line':
            self.view.setStyle({'model': self.proteinModelId},{'line':{coloring: color}});
            
        elif style == 'cartoon':
            self.view.setStyle({'model': self.proteinModelId},{'cartoon':{coloring: color, 
                                                                          'arrows': 'true'}
                                                              })
        else:
            self.view.setStyle({'model': self.proteinModelId},{'cartoon':{coloring: color,
                                                                          'arrows': 'true', 
                                                                          'tubes' : 'true'}
                                                              })
            
    def AddLigandWithStyle(self):
        """ add ligand and energy minimized ligand in viewer (called in render3D function) """
        if self.molViewState.ligandDict is not None:
            self.modelId = -1
            for models in self.molViewState.selectedModels:
                for modelName in list(models.keys()):
                    if modelName in self.childKeysForMolRender:
                        wgName = 'wg_'+ modelName + 'Visible'
                        if wgName in self.__dict__:
                            if self.__dict__[wgName].value:
                                self.view.addModel(models[modelName], 'sdf')
                                self.modelId = self.modelId + 1
                                self.AddLigandStyle(modelName, self.modelId)
                        else:
                            self.view.addModel(models[modelName], 'sdf')
                            self.modelId = self.modelId + 1
                            self.AddLigandStyle(modelName, self.modelId)
                            
    def AddProteinWithStyle(self):
        """ add protein in viewer (called in render3D function) """
        if self.molViewState.protein is not None:
            if 'wg_proteinVisible' in self.__dict__:
                if self.wg_proteinVisible.value:
                    pdb = Chem.MolToPDBBlock(self.molViewState.protein)
                    self.view.addModel(pdb,'pdb')
                    self.proteinModelId = self.modelId + 1
                    self.AddProteinStyle()
            else:
                pdb = Chem.MolToPDBBlock(self.molViewState.protein)
                self.view.addModel(pdb,'pdb')
                self.proteinModelId = self.modelId + 1
                self.AddProteinStyle()
                
    def render3D(self):
        """ This function updates the 3DMol.js viewer"""
        self.view.removeAllLabels()
        self.view.removeAllModels()
        self.view.removeAllSurfaces()
        
        self.GetWidgetValue()
        
        self.view.setBackgroundColor(self.wg_background.value)
        
        self.AddLigandWithStyle()
        
        # add label if required
        if self.labelPanel:
            self.AddLigandLabels()
            
        self.AddProteinWithStyle()
        
        # zoomTo does not work well for surface and label... so, zoomTo should not be default settings
        if self.onStart:
            self.view.zoomTo()
            self.onStart = False
            
        if self.propertyPanel:
            self.ShowLigandProperty()
            
        if self.emPanel:
            self.ShowMinimizationEnergy()
            
        display(self.view.update())
        
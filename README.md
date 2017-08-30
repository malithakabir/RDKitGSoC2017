# RDKitGSoC2017
## RDKit - 3Dmol.js integration
#### Mentors: Paul Czodrowski and Greg Landrum<br>Acknowledgement: Peter Gedeck reviewed codes, provided advice on code restructuring, and wrote initial MolViewState class. <br>Date: 29th August 2017<br>Email: malitha12345@gmail.com
#### To visitor: We appriciate your ideas on adding new features in it. So, please don't hesitate to drop your words. 
### Tasks listed in GSoC 2017
- A simple conformer browser (completed)
- Energy minimization of ligand extracted from PBD file (completed)
- Molecule editing (incomplete at the time of writing this - date : 29th August 2017)
### About Last 3 Months
#### Accomplishment
- We have created a maintainable codesbase. (links at the bottom of this page)
- We have developed 3D molecule viewer (that includes conformer (ligand) browsing) for RDKit.
- We have utilized RDKit "standard Mol object" (created by reading PDB file) to separate protein and ligand.
- We have enabled the viewer to show properties from SD file and minimization energy through ipython widgets.
#### Difficulties
- Persisting the 3D view angle while interacting with widgets
- Defining a data structure
- Making the long codes readable and maintainable
#### Notes On Pending Task
- We started getting output from the begining. The first two months we accomplished the tasks as written in proposal. And then we decided to improve code quality to make it maintainable. During the final month we have rewritten the codes several times with slight improvement of features. Therefore, time did not allowed us to accomplish 3rd task of the proposal.
- As I (Malitha Humayun Kabir) will keep working after the GSoC time frame, I will seek time from mentors regarding the 3rd task. So, we hope to get the 3rd task completed soon.
### Project Summary
- We are providing a function named "Check3Dmolpromise" that checkes whether the requied JavaScript libray "3DMol.js" loaded.
- We are providing a function named "ProcessLigandDict" that adds additional dict keys in a dictionary of molecule to allow 3D visualization.
- Our visualization system is a completely isolated python class named "MolView3D" that got the codes required for widgets rendering, widgets update, and molecule rendering persisting most recent visualization state.
- To feed molecule related data in "MolView3D", we have developed another completely isolated class named "MolViewState" that controls/decied what data will be rendered.
- The classes - "MolView3D" and "MolViewState" are NOT responsible for any kinds of calculation. Therefore, the property calculation and energy minimization system has been kept separate from those classes. Two functions - "MinimizeLigand" and "AddPropToLigandDict" does the energy minimization and porperty calculation respectively. 
- If viewer has been instantiated by view = ipy.MolView3D( necessary arguments here ) where ipy is the script named IpythonConsoleIntegration.py then the selected ids of the current visualization state can be accessed through view.molViewState.idPaired which is a tuple where each element of the tuple follows a systax (molecule Id, conformer Id). Note: since ids of the molecule and conformers can be accessed, user can export molecules through writing scripts (a notebook showing how to write molecules will be uploaded soon).
- We are also providing a function named "ExtractMolFragment" that can extract ligand from RDKit mol object. The mol object is expected to be created by reading PDB file using RDKit function named "MolFromPDBFile". User need to provide the residue name of the ligand to extact that from mol object using "ExtractMolFragment".
### Notes On MolView3D
#### Input Data
- Any RDKit mol object can be viewed after putting it inside a dictionary having molecule ids as keys of dictionary. Then "ProcessLigandDict" function will process the dictionary. The processed result then can work as input for "MolView3D".
- Except the processed dictionary of ligand and mol object of protein, there are only a few arguments that can be utilized by advanced users. 
#### Molecule Representations And Colors
- Representation for protein : cartoon, surface, line
- Representation for ligand (parent): line, cross, stick, ball-stick, surface
- Representation for ligand (energy minimized): line, cross, stick, ball-stick, surface
- Ligand color for both 'parent' and 'energy minimized' can be changed.
#### Mandatory Items And Panels In Viewer
- Molecule and conformer selection panel is mandatory for viewer. This panel can be minimal or full (inclued additional selection scheme).
- Visibility panel (checkbox name "ligandVisible", "proteinVisible", "emLigandVisible" (em = energy minimized)) is also mandatory.
- Since molecule viewer background color has impact on visualization clarity it is also madatory item.
- A button named "zoomTo" that used to focus the whole molecule has been kept as mandatory item.
#### Optional Panels In Viewer
- stylePanel - Panel for molecule representation
- propertyPanel - Panel for showing the properties exists in SD file and also the newly calculated properties
- emPanel - Panel for showing minimization energy if exists in supplied dictionary
- labelPanel - Panel for labeling conformers (molecule id and conformer id) and atoms (atom symbol and atom id)
### Submitted To Mentors For Evaluation
#### Through GitHub pull request to merge at rdkit:GSoC2017
- Notebook demonstrating usability : https://github.com/malithakabir/RDKitGSoC2017/blob/master/BrowseMultimolsV8.ipynb
- Visualization codes (Acknowledging Peter Gedeck): https://github.com/malithakabir/RDKitGSoC2017/blob/master/IPythonConsoleIntegration.py
- Ligand extraction codes (mainly written by Greg Landrum): https://github.com/malithakabir/RDKitGSoC2017/blob/master/LigandExtract.py
#### Not in pull request
- Demonstrating panels and conformer browsing : https://github.com/malithakabir/RDKitGSoC2017/blob/master/GSoC2017_notebook_1_ConformerBrowse_panels_and_confSelection.ipynb
- Demonstrating protein - ligand rendering and extraction of ligand from rdkit mol object : https://github.com/malithakabir/RDKitGSoC2017/blob/master/GSoC2017_notebook_2_ConformerBrowse_with_proteins.ipynb
- Demonstrating energy minimization and dictionary key utilization : http://34.210.39.113:8888/notebooks/GSoC2017_notebook_3_energy_minimization_overview.ipynb
- Demonstrating energy minimized ligand rendering with protein and not energy minimized ligand : http://34.210.39.113:8888/notebooks/GSoC2017_notebook_4_energy_minimization_and_visualization_(protein-ligand).ipynb
- A powerpoint slide (draft version) for RDKit UserGroupMeeting (will be held on September 2017).
### GitHub Pull Request Link
- https://github.com/rdkit/rdkit/pull/1484
- As per decision by mentors, codes from the aforementioned pulled request (all the development are in single pull request) would be merged today afternoon (european time) ( 29th August 2017)
### Contributors
- Paul Czodrowski was involved in answering questions, writing proposal, defining directions, communicating with Peter Gedeck and Brian Kelly for code review and merge respectively, and final work product submission. He offered a lot of disscussion over google handout, emails, and gmail's chat service. 
- Greg Landrum was involed in deciding how the code should be written and he appeared at every critical decision/time. He wrote the ligand extarction system and asked to me to modify that. 
- Peter Gedeck reviewed codes several times, wrote inital implementation of MolStateView class. His directions has made the codes maintainable.
- David Koes solved py3Dmol related issues with detail example.
- Brian Kelly reviewed for code mergeability.
- Malitha Humayun Kabir wrote most of the codes and was the participant of Google Summer of Codes (GSoc) 2017.
### Remark
- This page might not be updated in future for reference purpose but the repository may be updated with new notebooks (tutorials). 
### We appriciate your time here. Thank you very much for looking in it. Have a great day!

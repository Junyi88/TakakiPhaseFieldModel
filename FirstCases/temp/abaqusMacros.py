# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

def Macro1():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    session.mdbData.summary()
    o1 = session.openOdb(
        name='Y:/GitKrakenRes/TakakiPhaseFieldModel/FirstCases/G49/Polycrystal.odb')
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
        CONTOURS_ON_DEF, ))
    odb = session.odbs['Y:/GitKrakenRes/TakakiPhaseFieldModel/FirstCases/G49/Polycrystal.odb']
    nf = NumberFormat(numDigits=8, precision=0, format=ENGINEERING)
    session.fieldReportOptions.setValues(numberFormat=nf)
    session.writeFieldReport(fileName='SDVNodal.rpt', append=ON, 
        sortItem='Node Label', odb=odb, step=0, frame=87, outputPosition=NODAL, 
        variable=(('SDV26', INTEGRATION_POINT), ('SDV54', INTEGRATION_POINT), (
        'SDV56', INTEGRATION_POINT), ))
    odb = session.odbs['Y:/GitKrakenRes/TakakiPhaseFieldModel/FirstCases/G49/Polycrystal.odb']
    session.writeFieldReport(fileName='SDVIntegrationPoint.rpt', append=ON, 
        sortItem='Element Label', odb=odb, step=0, frame=87, 
        outputPosition=INTEGRATION_POINT, variable=(('SDV26', 
        INTEGRATION_POINT), ('SDV54', INTEGRATION_POINT), ('SDV56', 
        INTEGRATION_POINT), ))
    odb = session.odbs['Y:/GitKrakenRes/TakakiPhaseFieldModel/FirstCases/G49/Polycrystal.odb']
    session.writeFieldReport(fileName='Displacements.rpt', append=ON, 
        sortItem='Node Label', odb=odb, step=0, frame=87, outputPosition=NODAL, 
        variable=(('U', NODAL, ((COMPONENT, 'U1'), (COMPONENT, 'U2'), (
        COMPONENT, 'U3'), )), ))



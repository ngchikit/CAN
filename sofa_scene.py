

import Sofa.Core
import Sofa.constants.Key as Key 
from stlib3.physics.deformable import ElasticMaterialObject
import numpy as np

from scipy.optimize import least_squares

import Sofa
import time


class LinearForceAlternator(Sofa.Core.Controller):
    def __init__(self, node, points, force):
        super().__init__()
        self.node = node
        self.points = points
        self.force = force
        
        # Initialize force and times
        self.force_vectors = [[1, 0, 0]]  # Initial force as a list of numbers
        self.update_times = [0.3]  # Start with the first update time
        
        # Create the LinearForceField
        self.force_field = node.addObject(
            'LinearForceField', 
            points=points, 
            forces=[1, 0, 0], 
            force=force, 
            times=[0.3]
        )
        
        self.direction = 1
        self.time_interval = 0.3  # Time interval for force update
        self.last_update_time = self.node.getRoot().time.value

    def onAnimateBeginEvent(self, event):
        current_time = self.node.getRoot().time.value
        
        # Check if it's time to update
        if current_time - self.last_update_time >= self.time_interval:
            # Alternate force direction and append new force
            if self.direction == 1:
                self.force_vectors.append([-1, 0, 0])  # Append new force as a list
            else:
                self.force_vectors.append([1, 0, 0])  # Append new force as a list
            
            self.direction *= -1  # Toggle direction
            self.last_update_time = current_time

            # Calculate and append the next update time
            next_time = current_time + self.time_interval
            self.update_times.append(next_time)

            # Update the force field properties with list data
            self.force_field.forces = self.force_vectors
            self.force_field.times = self.update_times

def print_positions(root, indices):
    # Retrieve the positions of endo's degrees of freedom at specified indices
    endo_positions = np.array([root.ElasticMaterialObject2.dofs.position.value[idx] for idx in indices])
    print(f"Positions of endo.dofs at indices {indices}: {endo_positions}")
    return endo_positions

def fit_circle_2d(points):
    # Compute the mean of the points in the yz-plane
    y_m, z_m = np.mean(points, axis=0)

    # Calculate the circle parameters
    def calc_radius(c):
        return np.sqrt((points[:, 0] - c[0])**2 + (points[:, 1] - c[1])**2)

    def objective(c):
        return calc_radius(c) - np.mean(calc_radius(c))

    # Estimate the center using least squares
    from scipy.optimize import least_squares
    result = least_squares(objective, x0=[y_m, z_m])

    center_y, center_z = result.x
    radius = np.mean(calc_radius(result.x))

    return center_y, center_z, radius

# Define indices of interest
indices_of_interest = [0, 2, 4, 6, 9, 323, 325, 452]

def createScene(rootNode):
    from stlib3.scene import MainHeader
    MainHeader(rootNode, plugins=["SoftRobots", "Sofa.GL.Component.Rendering3D", "Sofa.Component.ODESolver.Backward",
                              "Sofa.Component.LinearSolver.Direct", "Sofa.Component.IO.Mesh", "Sofa.Component.Mass",
                              "Sofa.Component.SolidMechanics.FEM.Elastic", "Sofa.Component.Constraint.Lagrangian.Correction",
                              "Sofa.Component.MechanicalLoad"]
                              ,gravity=[0, 0, 0])
    rootNode.VisualStyle.displayFlags = "showVisual" # for rendering
    
    ################################ contact setting ########################################
    from stlib3.scene import ContactHeader
    rootNode.addObject('CollisionPipeline')
    ContactHeader(rootNode, alarmDistance=2, contactDistance=0.5, frictionCoef=0.1)
    ################################ contact setting ########################################

    MainHeader(rootNode, gravity=[0 ,0, 0],dt=0.01)

    sto = rootNode.addChild(ElasticMaterialObject(name='ElasticMaterialObject1', volumeMeshFileName="mesh/sto4_tran.vtk",
                                        translation=[0, 0, 0], scale=[8,8,8],
                                        surfaceMeshFileName="mesh/sto4_tran.obj",
                                        collisionMesh="mesh/sto4_tran.obj",
                                        surfaceColor=[1.0, 0.0, 0.0, 0.2],
                                        rotation = [0,0,0], youngModulus = 3000, poissonRatio = 0.3, totalMass = 10))

    sto.addObject('FixedLagrangianConstraint', name="FixedConstraint", indices="85 72 107")

    indices = set(range(652))
    # Remove the unwanted indices
    unwanted_indices = {85 ,72 ,107}
    indices.difference_update(unwanted_indices)
    sto.addObject('Monitor',
                       name='moni',
                       template='Vec3d',
                       listening=True,
                       indices=list(indices),
                       showPositions=False,
                       PositionsColor=[1, 1, 0, 1],  # Yellow color with alpha
                       showVelocities=False,
                       VelocitiesColor=[1, 1, 0, 1],  # Yellow color with alpha
                       ForcesColor=[1, 1, 0, 1],      # Yellow color with alpha
                       showMinThreshold=0.01,
                       TrajectoriesPrecision=0.1,
                       TrajectoriesColor=[1, 1, 0, 1],  # Yellow color with alpha
                       sizeFactor=0.1)
    
    sto.addObject('Monitor',
                    name='moni2',
                    template='Vec3d',
                    listening=True,
                    # indices=list(indices),
                #    indices=[267, 84],
                    indices=        [
                        474, 469, 468, 410, 409, 380, 381, 370, 371,
                        477, 408, 478, 407, 386, 383, 382, 384,
                        406, 405, 402, 385, 387, 391, 639
                    ],
                    showPositions=True,
                    PositionsColor=[1, 1, 0, 1],  # Yellow color with alpha
                    showVelocities=False,
                    VelocitiesColor=[1, 1, 0, 1],  # Yellow color with alpha
                    ForcesColor=[1, 1, 0, 1],      # Yellow color with alpha
                    showMinThreshold=0.01,
                    TrajectoriesPrecision=0.1,
                    TrajectoriesColor=[1, 1, 0, 1],  # Yellow color with alpha
                    sizeFactor=10)
    
    sto.addObject(LinearForceAlternator(node=sto, points="267 84", force="1000"))     # Attach the custom force alternator controller
    ####################### force for dynamics ###########################
    
    endo = rootNode.addChild(ElasticMaterialObject(name='ElasticMaterialObject2', volumeMeshFileName="mesh/endo1.vtk",
                                    surfaceMeshFileName="mesh/endo1.obj",
                                    collisionMesh="mesh/endo1.obj",
                                    surfaceColor=[0.0, 1.0, 0.0, 1],
                                    translation=[32.0, 0.0, 0.0], rotation = [0,0,0], youngModulus = 15000, poissonRatio = 0.3, totalMass = 10, scale = [8,8,8]))
    

    ########################## Add fix box at the end of endo #############################
    translatebox = [12,0,0] # it should be changed by considering the translation of endo
    atPositions=[-33,-10,-10,-23,10,10]
    atPositions[0] += translatebox[0]
    atPositions[3] += translatebox[0]
    from stlib3.physics.constraints import FixedBox
    FixedBox(endo, doVisualization=True, atPositions=atPositions)


    ########################### todo: add monitor to visulize force ###################
    # Create the set of indices
    indices = set(range(559))

    # Remove the unwanted indices
    unwanted_indices = {133, 258, 304, 389}
    indices.difference_update(unwanted_indices)

    ################################# initialize step controller ###############################
    # TODO: initialize step controller
    from controller import EndostepController
    endo.addObject(EndostepController(endo))

    ################################# initialize step controller ###############################

    ################################## see a point in the scene ############################
    # position = np.array(endo.dofs.position.value[390])
    tiptranslate = [12.0,0,0]
    tipposition = [30.55632019, 0.4443016052246096, -0.2793472707271626]
    tipposition[0] += tiptranslate[0]
    tip1=endo.addChild('Tip1')
    tip1.addObject('MechanicalObject',name='tip1', position=tipposition
                   ,showObject = True, showObjectScale = 50) # position after adjusting
    tip1.addObject('BarycentricMapping', mapForces=False, mapMasses=False)
    
    from softrobots.actuators import PullingCable
    from splib3.loaders import loadPointListFromFile
    from cablecontroller import cablestepcontroller
    from cablecontroller import cablekeyboardcontroller
    ######################################## Todo cable constraints and cable control ########################
    cabletranslate = [12,0,0] # it should be changed by considering the translation of endo
    cable1 = PullingCable(endo,
                        "cable1",
                        # pullPointLocation=[18.55632019, 1.4443016052246096, 0.7206527292728374],
                        # rotation=rotation,
                        translation=cabletranslate,
                        cableGeometry=loadPointListFromFile("mesh/cabledesign/cable1.json"));
    cable2 = PullingCable(endo,
                        "cable2",
                        # pullPointLocation=[18.55632019, 1.4443016052246096, -1.2793472707271625],
                        # rotation=rotation,
                        translation=cabletranslate,
                        cableGeometry=loadPointListFromFile("mesh/cabledesign/cable2.json"));
    cable3 = PullingCable(endo,
                        "cable3",
                        # pullPointLocation=[18.55632019, -0.5556983947753904, 0.7206527292728374],
                        # rotation=rotation,
                        translation=cabletranslate,
                        cableGeometry=loadPointListFromFile("mesh/cabledesign/cable3.json"));
    cable4 = PullingCable(endo,
                        "cable4",
                        # pullPointLocation=[18.55632019, -0.5556983947753904, -1.2793472707271625],
                        # rotation=rotation,
                        translation=cabletranslate,
                        cableGeometry=loadPointListFromFile("mesh/cabledesign/cable4.json"));
    endo.addObject(cablestepcontroller(cable1,cable2,cable3,cable4))
    ######################################## Todo cable constraints and cable control ########################

    #################################### monitor contact force #########################################
    # Todo: Contact force monitoring

    listener = rootNode.addObject(
        "ContactListener",
        collisionModel1=sto.CollisionModel.LineCollisionModel.getLinkPath(),
        collisionModel2=endo.CollisionModel.LineCollisionModel.getLinkPath(),
    )
    #################################### monitor contact force #########################################
    from splib3.animation import AnimationManager

    return rootNode
USE_GUI = True

def main():
    import SofaRuntime
    import Sofa.Gui
    # Make sure to load all SOFA libraries

    #Create the root node
    root = Sofa.Core.Node("root")
    # Call the below 'createScene' function to create the scene graph
    createScene(root)
    Sofa.Simulation.init(root)


    ####################################### Cable todo ######################################

    # Print the positions of specified indices
    positions = print_positions(root, indices_of_interest)

    # Extract y and z coordinates for circle fitting in the yz-plane
    yz_points = np.array([(y, z) for _, y, z in positions])
    # Perform the fit and print the results
    center_y, center_z, radius = fit_circle_2d(yz_points)
    print(f"Center in yz-plane: ({center_y}, {center_z}), Radius: {radius}")

    ####################################### Cable todo ######################################

    ############# Direct start animation after initialize the rootNode ###################
    # Todo2: start animation
    root.animate = True
    ############# Direct start animation after initialize the rootNode ###################

    ################################ Get the data ########################################
    # Todo3: get the data
    # State-based DRL need the info of EE & target & force
    def state(root):
        target = np.array(root.point_tar.point_tar.rest_position.value)  # Assuming rest_position contains the positions of all nodes
        print(f"Tar: {target}")
        EE = np.array(root.ElasticMaterialObject2.Tip1.tip1.position.value)
        print(f"Tip: {EE}")
        print("Contact indices", root.ContactListener.getContactData().get('collisionElementsModel2'))
        print("Contact force", root.ContactListener.getContactData().get('collisionPointsModel2'))
        print("Contact source", root.ContactListener.getInstanciationSourceFilePos())
    ################################ Get the data ########################################
    
    if not USE_GUI:
        # Real-time visualization and recording loop
        for iteration in range(512):
            print(f"In iteration {iteration}")

            # Apply action and animate the simulation
            root.ElasticMaterialObject2.EndostepController.applyAction([0.1, 0, 0])
            Sofa.Simulation.animate(root, root.dt.value)
            state(root)
            # Retrieve contact data
            contact_indices = root.ContactListener.getContactData().get('collisionElementsModel2')
            contact_forces = root.ContactListener.getContactData().get('collisionPointsModel2')

            # Update images
            update_images(contact_indices, contact_forces)

            # Resize images to fit display window
            grayscale_display = cv2.resize(grayscale_image, display_size, interpolation=cv2.INTER_NEAREST)
            rgb_display = cv2.resize(rgb_image, display_size, interpolation=cv2.INTER_NEAREST)

            # Display the images in configured windows
            cv2.imshow('Force Magnitude (Grayscale)', grayscale_display)
            cv2.imshow('Force Components (RGB)', rgb_display)

            # Write frames to video files
            grayscale_video_writer.write(grayscale_display)
            rgb_video_writer.write(rgb_display)

            # Refresh display and add delay for visualization
            cv2.waitKey(1)  # 1 ms delay to allow updating
    else:
        # Find out the supported GUIs
        print ("Supported GUIs are: " + Sofa.Gui.GUIManager.ListSupportedGUI(","))
        # Launch the GUI (qt or qglviewer)
        Sofa.Gui.GUIManager.Init("myscene", "qglviewer")
        # Sofa.Gui.GUIManager.Init("myscene", "qt")

        Sofa.Gui.GUIManager.createGUI(root, __file__)
        # 2560x1600
        Sofa.Gui.GUIManager.SetDimension(2560, 1600)

        # Initialization of the scene will be done here
        Sofa.Gui.GUIManager.MainLoop(root)
        Sofa.Gui.GUIManager.closeGUI()
        print("GUI was closed")

    print("Simulation is done.")
    # Release video writers to save the video files
    grayscale_video_writer.release()
    rgb_video_writer.release()
    # Close the display windows after the loop
    cv2.destroyAllWindows()
    print("Videos saved successfully.")
# Function used only if this script is called from a python environment
if __name__ == '__main__':
    main()

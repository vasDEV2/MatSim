### READ THIS BEFORE USE

This setup will first open a map of waypoints from which waypoints have to be picked; This will then
open up the simulator and run the simulation. Once the simulation is done, matlab auto-exits. VDM is
consequently initialized producing a comparison between the simulator and the VDM prediction. Once the 
whole process is finished, the generated data is automatically moved into a folder called DataSetX where
X is the latest index.

### STEPS TO FOLLOW:
 
1. Once the zip files are extracted, execute the "executesim.sh" file by typing "./executesim.sh".

2. Wait for the map to show up. once up, select the desired path and click on 'export to workspace". Manually 
   close the map to continue the script.
   
3. The simulator will start next, where the car will follow the selcted waypoints. You can finish the simulation 
   by either manually closing the simulation window or waiting for 300secs(compulsory to initiate next steps).
   
4. Once the simulation exits, the vdm scripts will run, plotting the map. after analyzing the map, you will need
   to close it for the bash script to proceed.
   
5. The script will finally create the "datasetx" folder into which the newly generated data is backed up. This
   will ensure that you don't need to backup the files before running the next simualtion.
   
### KNOWN ISSUES  
 
 1. I have tried to make all the paths relative but incase I have missed any, let me know.
 
 2. High frequency noise in the steering signal.
 
 3. The waypoint map is quite inaccurate and you need to do some guessing around to make it follow the correct path.
    highly likely that the car will roughly follow the selected path.
    
 4. Matlab was refusing to quit in the bash script so I have created this weird flag system for it to exit. Let me know
    if this breaks.


0. If you haven't compiled the code in your system, compile it by:

chmod u+x mk.sh

./mk.sh

This step can be skipped once code is compiled.

1. Inside the PntWrks directory create a folder called "sim"

2. Copy a template test case from the "templates" folder and rename folder to "in" and place it inside the "sim" directory

3. edit the input parameters for the simulation as needed and run the simulation by:

chmod u+x run.sh 

./run.sh

4. To use the 2D python plotter:

python3 plotpnts.py <variable 1> <variable 2> <scale val> <point size>

for example:

python3 pltpnts.py u phi 5 3

Or to use the 2D meshed data plotter:

python3 pltmsh.py <variable 1> <variable 2> <scale val> <show mesh flag>

for example:

python3 pltmsh.py u phi 5 0

3D data is best viewed by opening the generated vtk files in paraview or VisIt.

5. To export the generated simulation pics to a video run script:

chmod u+x genvid.sh

./genvid.sh <variable name>

for example:
./genvid.sh phi

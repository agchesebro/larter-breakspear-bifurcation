# Running MatCont analyses from the manuscript

To run the analyses from Section 3 and generate Figures 2a, 3a, and 4a, follow the steps below.

1. Download and install MatCont following the instructions from the MatCont developers.
2. In the main MatCont folder, find the subfolder Systems. Place the contents of the Systems folder from this repo into that folder as well.
3. Open Matlab and run MatCont
4. Choose Select > System > Load/Edit/Delete Systems.
5. Select larter-breakspear from the list on the left and click the Load button.
6. In the Main window, select Type > Initial Point > Point.
7. This should load the variables into the Integrator and Starter windows. Check to make sure the parameters all match Table 1 (and t=0, V and Z are in (-1, 1), and W is in (0, 1)). 
8. To observe the parameters changing, in the main window select Window/Output > Numeric (launches Numeric window) and Window/Output > Graphic > 3D plot.
9. Select Compute > Forward. This will produce a graph of the standard oscillations in 3D.
10. Once this has run, select Type > Initial Point > Equilibrium.
11. In the Starter window, you need to select the ion you’re interested in continuing (V_Na, V_K, or V_Ca).
12. In the Main window, select Compute > Forward (may need to do Compute > Backward as well depending on your initial conditions).
13. Find the Hopf point (the Control window will pause, and you can click Stop).
14. Double click on the Hopf point in the Graphic window and say Yes to the prompt to load it.
15. In the Main window, select Compute > Forward and you’ll generate the results for the figure. Whenever it pauses at a point, note the values in the Numeric window.

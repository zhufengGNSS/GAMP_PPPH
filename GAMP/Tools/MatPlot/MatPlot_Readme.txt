1. The file named "analysis.ana" will be generated after pushing the button of "Analysis" on the GUI of MatPlot.
The format of "analysis.ana" is:
the 1st column is site name;
the 2nd column is convergence time (min) in east component;
the 3rd column is convergence time (min) in north component;
the 4th column is convergence time (min) in up component;
the 5th column is convergence time (min) of PPP processing for the site;
the 6th column is positioning accuracy (cm) in east component before convergence;
the 7th column is positioning accuracy (cm) in north component before convergence;
the 8th column is positioning accuracy (cm) in up component before convergence;
the 9th column is positioning accuracy (cm) in east component after convergence;
the 10th column is positioning accuracy (cm) in north component after convergence;
the 11th column is positioning accuracy (cm) in up component after convergence.

2. The label of ¡°ylim(ppp)¡± on the GUI of MatPlot can be used to set Y-axis range. 
Note that ¡°ylim(ppp)¡± sets the maximum value along the axis and the negative of this value is the minimum along the axis. 
A value of 10, for example, plots the Y-axis as -10 to 10. Setting 0 uses a MATLAB default along the Y-axis. 
It is recommended to set this value greater than or equal to 0.

3. The executable version of MatPlot is MatPlot.exe.
But the ¡°MCR_R2016b_win64_installer.exe¡± should be firstly download from
https://www.mathworks.com/products/compiler/matlab-runtime.html .
Doulbe click and install it before you run MatPlot.exe.


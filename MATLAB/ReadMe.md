Before running anything, it is important to import all the necessary files.
This can be done by running the config.m file.

In the folder tests, there are multiple scripts that can serve as an example on how to run the code.
For instance we can run the script test_disc_rho_models.m that is modelling a discontinuous initial density profile.

Once we run a test script, we can compute the relative error by running the compute_error.m script (for tests comparing collision rates) of compute_error_mom.m (for tests comparing multiple moments). NOTE: It is necessary to change both scripts depending on the test case. 

Once the errors are computed, these can be plotted using the script plotter_error_models.m or plotter_error_mom.m in the folder "plotter". Similarly, we can also plot the computed solutions using the scripts plotter_models.m or plotter_moments.m

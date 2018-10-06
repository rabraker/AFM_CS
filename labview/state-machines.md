
## Notes about different state machines for different applications.
When doing raster-scans, pure CS (no mu-paths) or CS-mu, there are different states which need to be implemented as state machines. (These states are not to be confused with the states of a dynamical system...)

* Raster scan states.

 1. Lift tip after disengaging agilent PI control
 * Lower tip. Deflection signal remains constant until tip engages the surface. Detect a change greater than some threshold to trigger the next state.
 * PID control. State reading raster trajectory from FIFO buffer. Engage PID in z-axis, using the control input from state 1 as the initial condition for the integrator. Write data to fpga-to-host FIFO. Trigger to state 4 if (a) either fifo buffer has timed-out, (b) a boolean value (indicator PID-to-tip-up) is set to true.
 * Currently, exit the main control loop and initiat graceful return to the origin.
 
* Pure CS-states
  This description assumes we always measure u_x=u_y=0. This is of course not necessary, but I think it simplifies the logic.
    1. Initial Descent. Effectively the same as raster-scan descent I think. Deflection signal remains constant until tip engages the surface. Detect a change greater than some threshold to trigger the next state.
    * PID-control on z-ax (and possibly also on x&y). Measure z-axis height. 
    * Tip up. Our AFM has no sensor, so the only thing we can do is give a step command and wait a certain number of clock ticks. Yufan may be able to do something a little more interesting. 
    * x-y-axis move. Trigger on detection of settling, or use heuristic, or ???. Transition to state 1.


* $\mu$ path CS-states

    1. Initial Descent. Effectively the same as raster-scan descent I think. Deflection signal remains constant until tip engages the surface. Detect a change greater than some threshold to trigger the next state.
    * Initiate PID-control on z-ax. Transition to PID on x&y. Trigger on settling or clock-tick hueuristic.
    * Start mu-path scan. I.e., start reading trajectory data from FIFO and track. Send z, err, x, & y data back to host. Trigger on number of samples. (For only mu-paths, we probably don't need to send the trajectory to track from the host. As I understand it, this will be a simple ramp signal in one or both directions. If we hold the ramp rates constant, we only need to send a clock tick count (integer) from the host for each mu-path. On the other hand, if we did send the whole trajectory from the host, we would be able to follow arbitrary patterns, and the code for mu-paths would work the same for anything else (modulo the control law being the same). 
    * Tip up. Our AFM has no sensor, so the only thing we can do is give a step command and wait a certain number of clock ticks. Yufan may be able to do something a little more interesting. 
    * x-y-axis move. Trigger on detection of settling, or use heuristic, or ???. Transition to state 1.


From a hardware efficiency standpoint, it doesn't make sense to repeat functionality inside each state of the state machine. For instance, if a PID controller is required in more than once state and we repeat that logic inside multiple cases of the case structure, then we (a) increase hardware usage, (b) make bugs more likely and debugging more difficult.  Thus, I think it makes sense to limit the logic inside inside the state machine only to decision making. Other logic/control should be kept outside the main state machine and it's computations either disgarged, or its execution enabled/disabled with a boolean controlled by the state machine. 

# Data logging.
In the last two schemes, it's not obvious to me which states should be logging data. The easiest thing to do is to send all 4 measurements back all the time. But then I think post processing is going to be a nightmare, figuring out which data come from which state. For instance, the data that comes from the 'xy-move' state is pretty irrelevent as far as I know. So maybe the thing to do is pack another value, an 8bit integer along with the data which tells which state it came from. 


## the data indexing problem
The reason this is so difficult is that if we acquire data and sent it back to the host ALL THE TIME, some of that data is meaningless for image reconstruction. Ie, when we move from point-to-point, or mu-scan-to-mu-scan, the data in those transitions doesn't help get an image. 

Essentially, we need to send meta data with the data. 

This is because between the fpga and RT or host computer, we don't have any kind of real time communication mechanism. AFAIK, they are all lossy, with the exception of the FIFO, and none will occure in RT.


I propose to solve this problem by what I will refer to as *measurement entities* (MI). An MI might be a single pixel measurement point, a mu-scan etc. But in general, we will have N MIs. For every CS modality, the respective MIs chare common characteristics:
1. x&y start coordinates
* an index, ie, MI number 235
* some kind of metric of how "long" we should spend at that MI.

The crucial point is, that if we package the MI-index with the x,y,z,err, u data in the fifo stream, post processing becomes a LOT more straightforward. 

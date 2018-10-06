# Basic Implementation of general CS Scanning

The purpose of this document is to provide a description and basic usage help for my implementation of CS in labview, which is contained in the labview project `afm_imaging_control.lvproj`. The main host file to run is `play-AFM-CS-imaging.vi`. 

Running this vi requires a number of front panel controls to be set as well as a trajectory `.csv` file. Before describing what all those controls do, I first need to describe the underlying fpga vi that the host vi controls, which is `cs-scanning.vi`. 

`cs-scanning.vi` implements a single control loop, i.e., with one sample frequency for all three axes. I am using a sample rate of 25khz, which at a 40mhz FPGA clock rate is 1600 clock ticks. The behaivior of the control is governed by what is essentially a state machine, implemented as a set of case structures. In general, the PID control runs continusously on the x and y axes. The setpoints or trajectories to follow are determined by the state machine. The z-axis transitions between closed loop and open loop control and these transitions are governed by the state machine.

These state are all triggered sequentially, except for the transition from 5 back to 1. They are:

0. User initialization. You essentially always want to start in this state. The xy setpoint is nominally 0 (more on this in a bit). For our setup, you would turn off the agilent PI controller, and start up this vi. You then manually pull the tip away from the surface with the fron panel control `init Z up`. 
1. xy move and detect steady state. In the first iteration, we read an xy setpoint from the FIFO buffer and set x_ref and y_ref to those values. In subsequent iterations we do NOT read from the FIFO buffer, but just recycle that same setpoint. The trigger to move to the next state is a detection that the x(k) and y(k) have reached a settling criterion.
2. Tip lower. We now lower the tip towards the surface, at a rate determined by rate `ramp-rate`. Once the absolute value of the deflection (or error) signal exceeds `TOL`, that is the trigger to move to the next state.
3. We now start PID control on the z-axis and wait for `|deflection(k) - setpoint|` to reach a settling criterion, which is the trigger to move to state 4.
4. Start measuring the surface. We initiate reading a trajectory path from the host-to-FPGA FIFO buffer, as well as start logging data into the fpga-to-host FIFO buffer. The xy-axis are now tracking some kind of trajectory, but our control here is agnostic to wheather this is a my-path, spiral scan, or discrete point. We keep on doing PID control on the z-axis. The trigger to move to the next state is that the xy trajectory we are following ends. I determine this by packing meta data into the host-to-FPGA FIFO data, which I will describe more fully later.
5. Tip up and wait. The xy PID control is set to the last value of the trajectory we were following, and we issue a step up command to the z-axis control. For my hardware, I just wait `zup-N` samples, since I can't measure anything once the surface is disengaged. After `zup-N` samples we transition back to state 1.

## FIFO data packing specifications
For all of this to work properly, data needs to packed into both `host-to-FPGA` and `fpga-to-host` in a specific format. One of the problems I was trying to solve with this is how do you determine what is a new setpoint and what is a trajectory? 
### host-to-FPGA
Both setpoint data (ie, move-entities) and measurement trajectory (ie, measurement entitities) come into the FPGA via the same fifo buffer. To distinguish which is which, I pack meta data with it, which is just a third number, which is always some integer. For new setpoint or movement data, this is always `j=0`. Then, when reading a new setpoint, we expect FIFO data to look like

```
[x_r(1), y_r(1), 0]
```

In my current setup, I expect *one* set that looks like that. The next set should correspond to a CS measurement. It will look like

```
[x_r(1), y_r(1), j, x_r(2), y(r), j, ... x_r(N), y_r(N), -1]
```

Notice that the index (meta data) j is the same, because j should be the index of the CS point we are currently measuring. The only exception to this is the end: when we reach the end of the current CS trajectory to follow, we expect `j=-1`, which is our trigger to stop measuring and transition from state 4 to state 5. 

If we timeout on reading the FIFO buffer and miss data, this whole scheme falls apart. That is why my vi is set to abort if any FIFO timout occurs.
### fpga-to-host
Currently, I am logging data for both states 1 and states 4, i.e., both moves and measurements. In reality, we could probably drop logging data for moves. The logging FIFO spec is similar, but we have more data to pass:

```
[x(1), y(1), e_z(1), z(1), u_z(1), j, x(2), y(2), e_z(2), z(2), u_z(2), j,...]

```

Again, we differentiate move from measurement data. For moves, we set j=0, while for measurements we set j equal to the CS index we got from `host-to-FPGA`. This, I believe, should ease post-processing on the host side. 

## Miscellanouse Concerns
### Instability detection
For all three axes, I check at each sample period if the control input exceeds some bound. If this is the case, I assume something has gone awry and abort. You could do something more sophiscticated I'm sure, but this works pretty well.

### FPGA exit
It is important to ensure the last values the FPGA writes to the DACs are all zero, because the DACs will hold onto their last value even after the FPGA vi exits. This is accomplished in the subvi `gracefull_return.vi`, which slowly returns all three controls towards zero, to avoid giving an uncontrolled large step input. I believe this vi should be hard coded, and not subject to user tweaking, which prevents accidently entering destructive values.

### XY nonzero initial condition.
On our piezo stage, the x and y sensor generally generate report a nonzero value for a zero control input. Therefore, I measure this initial offset before the main control loop starts, and subtract it off of all subsequent ADC measurements. 

### Functionality separation
I have tried to separate decision making from other logic in the state machine. For example, data logging, and FIFO buffer reading, and xy-axis PI control happen in more than one state. Thus, the code to do these tasks is moved outside the state machine, which prevents having code that does the same thing in more than one state. This has a couple of advantages:
1. You only have to update the code in one place.
2. This should decrease FPGA fabric consumption.

This could also be done for the z-axis PI control, but I haven't gotten around to it yet. 



## Host VI and Front Panel Controls

### Tab 1
These are the most commanly used settings.
1. `init-Zup`. Increase this to disengage the surface.
* `trigger-0`. Hitting this transitions us from state 0 to state 1.
* `z-UP`. This is the amount, in volts, to step the z-axis up in state 5. 
* `numEntities`. The total number of CS entitites we intend to visit. It would be better to read this value off of some meta data in a header of the data-in.csv file
* `ki-xy`. This is the control gain for the x and y PI controller transfer function, which looks like ki/(z -1). 
* `TsTicks` This sets the sample frequency. 1600 ticks at 40Mhz is 25khz.
* `TOL`. This is the tolerance to detect sample engagement in state 2.
* `setpoint` This is the z-axis setpoint, ie error-signal=(setpoint - deflection_signal)
* `z-u-max`. The maximum allowable z-axis input before we conclude instability and abort the entire thing. 
* `ramp-rate` How fast we approach the sample in state 2. This is in volts per sample period and MUST be negative (for our system). -5e-6 is really slow, -5e-5 is fast in human time. It would be good to make this as fast as possible. 

### Tab 2: xy-axis tweaks
1. `xy error threshold`. This determines how close to the setoint we must be to consider the x and y axes to have settled in state 1. 
* `xy-settled samples`. How many samples the x or y axis must be inside of the error threshold before we declare settling.

### Tab 3: z-axis tweaks
1. `z-error threshold`,This determines how close to the setoint we must be to consider the z axis to have settled in 3.
* `z-settled samples`. How many samples the z  axis must be inside of the error threshold before we declare settling. 
`z-up-N` How many samples do we wait for in state 5.


# Associated Matlab Classes
The input data is specified in a fairly simple format. However, I have written a set of matlab classes that help put this in a common framework for arbitrary measurement trajectories. More information on how to use them can be found in `matlab-code/readme.md` and `matlab-code/examples/cs-trajectories-example.m`.










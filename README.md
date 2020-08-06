# Radar target generation and detection

Final project repo from Radar module of Sensor Fusion Nanodegree (Udacity)

### Workflow :

1. Range 1D-FFT for range detection

2. Doppler 2D-FFT for velocity and range detection

3. CFAR clutter removal technique and thresholding

4. Final Project : apply all the above with 2-D CFAR
### Directory contents :

1. `Range_detection/` : FFT1D.m , range.m (range detection)

1-D FFT for range detection 

![alt-text](https://github.com/aditya-167/Radar-Target-Generation-Detection/blob/master/Output/1dfft.jpg)


2. `Doppler_velocity_detection` : FFT2D.m , dopplerVelocity.m

2-D FFT for doppler velocity and range detection plot

![alt-text](https://github.com/aditya-167/Radar-Target-Generation-Detection/blob/master/Output/2dfft.jpg)

3. `Clutter_and_thresholding` : CFAR1D.m (only 1-D CFAR here, 2-D CFAR in the project.m file)

1-D CFAR with self adjusting threshold based on noise and object density.

![alt-text](https://github.com/aditya-167/Radar-Target-Generation-Detection/blob/master/Output/CFAR1D.jpg)

##

### 2D CFAR process
* Loop over elements of RDM array each iteration selecting one cell to be the CUT (Cell Under Test)<br>
`for i = Tr+Gr+1 : (Nr/2)-(Gr+Tr)`<br>
`for j = Td+Gd+1 : Nd-(Gd+Td)`
* For each iteration loop over the training cells excluding the guarding cells to sum their values<br>
`for p = i-(Tr+Gr) : i+(Tr+Gr)`<br>
`for q = j-(Td+Gd) : j+(Td+Gd)`
* Calculate the average of the noise value<br>
`noise_level = noise_level + db2pow(RDM(p,q));`
* Convert using pow2db<br>
`th = pow2db(noise_level/(2*(Td+Gd+1)*2*(Tr+Gr+1)-(Gr*Gd)-1));`
* Add the offset value
* If the CUT is greater then the threshold replace it by `1`, else `0` 

### Training, Guard cells and offset
* `Tr = 12, Td = 10` For  Range and Doppler Training Cells.
* `Gr = 6, Gd = 6` For  Range and Doppler Guard Cells.
* `offset = 1.4` offset value.

### suppress the non-thresholded cells at the edges
sclicing the output such that we have the surrounding rows and columns depending on the Training cells for both range and doppler.<br>
`RDM(union(1:(Tr+Gr),end-(Tr+Gr-1):end),:) = 0;  % Rows`<br>
`RDM(:,union(1:(Td+Gd),end-(Td+Gd-1):end)) = 0;  % Columns`

### Final Results:
final output for a target at 50m moving at -25 m/s relative velocity<br><br>

Range detection 50m

![alt-text](https://github.com/aditya-167/Radar-Target-Generation-Detection/blob/master/Output/Range.jpg)

Range and velocity from 2D-FFT

![alt-text](https://github.com/aditya-167/Radar-Target-Generation-Detection/blob/master/Output/doppler.png)

2-D CFAR

![alt-text](https://github.com/aditya-167/Radar-Target-Generation-Detection/blob/master/Output/CFAR.jpg)

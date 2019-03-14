# Repeating_Ear
\chapter{Determination of the relative location between two subevents}

\section{Waveform cross-correlation to determine the precise differential arrival time}
We apply a waveform cross-correlation technique to estimate the precise arrival time difference between two phases. 
Waveform cross-correlation has been used extensively to investigate waveform similarity [e.g., Shearer, 1997; Waldhauser \& Ellsworth, 2000]. The similar focal mechanisms of two doublet earthquakes will result in highly similar waveforms. The normalized cross-correlation function between two signals x(t) and y(t) can be computed as: %We make the cross-correlation time-window one second longer than the actual durations of the signals and compute : 
\begin{align} \label{eq:1}
	\phi_{xy}^{N}(\tau) = \frac{\int^{+\infty}_{-\infty} x(t)y(t+\tau)\ \ dt}{\sqrt{\phi_{xx}(\tau)\phi_{yy}(\tau)}}
	\end{align}
where $\phi_{xx}$ and $\phi_{yy}$ denote the auto cross-correlation functions of x(t) and y(t), and $\tau$ denotes the time-lag. \\
If two signals are similar enough, the maximum amplitude of the cross-correlation function is well-defined. We can measure the time-shift $\tau$ that has the maximum cross-correlation coefficient and then use it to calculate the arrival time difference:
\begin{align} \label{eq:1}
	\Delta t = \Delta t_{0} + \tau, 
	\end{align}
where $\Delta t_{0}$ denotes the time difference between the onsets of two time windows.  As an example, Figure 3 shows the cross-correlation process at station YSS for the 1994 Kuril Earthquake. For this example, the maximum cross-correlation of 0.98 is for an arrival time difference of 1.95 seconds. The measurement errors of the differential arrival times can be smaller than the sampling rate [Waldhauser \& Ellsworth, 2000]. The sampling rates of our stations are either 20 samples per second or 40 samples per second. Therefore, the uncertainties of measurement are estimated to be $<$ 0.05 sec. 
\begin{figure}[hbt!]
\centering
\includegraphics[width=14cm]{tmp1.jpg} 
\caption{Waveform cross-correlation of a P-wave at station YSS. This waveform is from the 1994 Kuril Earthquake. a) Truncated data from one waveform with two time windows of similar phases indicated by the red dashed lines. b) The cross-correlation coefficient of two phases as a function of arrival time difference. c) Relative amplitude of two phases. The solid blue curve represents the time-series in the first time-window. The second red solid curve represents the time-series in the second time-window, which is shifted to the left by the differential arrival time got in the cross-correlation procedure. }
\label{Figure 1}
\end{figure}

\section{Forward-modeling of travel-time differences}
\begin{figure}[hbt!]
\centering
\includegraphics[width=14cm]{Figure_1.png}
\caption{Cross-section showing the geometry of the ray path from the earthquake to the station. The yellow stars denote the subevents and the red square denotes the seismic station. The green curve indicates the synthetic seismogram recorded by the seismic station if noise-free. The inset shows the source region of the two subevents. The subevents are separated by $\Delta x$, $\Delta y$, and $\Delta z$ km in the west-east, south-north, and vertical directions. Angle $i_0$ denotes the takeoff angle of the seismic ray. }
\label{Figure 1}
\end{figure}
We assume that the distance between two subevents is much smaller than the distance between the event and the station. Therefore, the only difference between the ray paths of two subevents is in the source region. Figure 3.1 shows the ray paths from the seismic event to the seismic station. We also assume that the velocity structure in the source region is homogeneous. Under these assumptions, the travel-time difference $\Delta t_{tra}$ between two subevents can be expressed as: 
 \begin{align} \label{eq:1}
	\Delta t_{tra} = \frac{-\Delta x \sin{i_{0}}\cos{az} - \Delta y \sin{az}\sin{i_0} + \Delta z \cos{i_{0}}}{v_p}, 
	\end{align}
where $\Delta x$, $\Delta y$, and $\Delta z$ denote the distance between two events in the west-east, south-north and vertical directions, $i_0$ denotes the takeoff angle of the seismic ray, az represents the azimuth from the source to the station and $v_p$ is the P-wave velocity at the depth of the source. The takeoff-angle of the seismic ray can be estimated as: 
 \begin{align} \label{eq:1}
	i_0 = \text{asin}(\frac{pv_p}{r}), 
\end{align}
where p denotes the ray parameter of the ray arriving at a given distance which can be estimated using the TauP software package [Crotwell \textit{et al.}, 1999] based on the IASP-91 model [Kennett and Engdahl, 1991], and r is the radius at the earthquake source. \\
\indent To show the corresponding pattern of predicted travel time difference between two subevents, we simulate four models using the P-wave velocity at the source depth of the 2000 Tonga earthquake (9.79 km/s). We set the spatial separation parameters $(\Delta x, \Delta y, \Delta z)$ to be (1 km, 0 km, 0 km), (0 km, 1 km, 0 km), (0 km, 0 km, 1 km), and (0.5774 km, 0.5774 km, 0.5774 km), respectively, so that for each model, the total spatial separation between the two subevents is equal to 1 km. The resulting patterns of the predicted travel time differences are shown in the Fig. 3.2.  In the model with the two events separated only in the east-wast direction, the maximum absolute time separation time will reached when the takeoff-angle is 90$^{\circ}$, and the azimuth is 0$^{\circ}$ or 180$^{\circ}$, and the minimum absolute time separation will be reached when the azimuth is 90$^{\circ}$ or 270$^{\circ}$. In the model with the two events separated only in the north-south direction, the maximum absolute time separation will be reached when the takeoff-angle is 90$^{\circ}$, and the azimuth is 90$^{\circ}$ or 270$^{\circ}$, and the minimum absolute time separation will be reached when the azimuth is 0$^{\circ}$ or 180$^{\circ}$. In the model with the two events separated only in the vertical direction, the maximum absolute time separation will be reached when the takeoff-angle is 0$^{\circ}$, and the minimum absolute time separation will be reached when the takeoff-angle is 90$^{\circ}$. In the last model, the maximum absolute time separation will be reached when the takeoff-angle is 60$^{\circ}$, and the azimuth is 220$^{\circ}$.

\begin{figure}[hbt!]
\begin{center}
\includegraphics[width=18cm]{Figure_2.png}
\caption{The patterns of the predicted travel time differences for four models with different spatial separations between subevents. Each row is labeled with the spatial separation. The first column of graphs show the spatial separation of the two subevents (red dots) with the first subevent at the origin. In the second column, the travel time differences are plotted on the lower-hemisphere stereonets. The colors in the stereonets represent the travel time differences. The third column of graphs indicate travel time differences as a function of station azimuths. The fourth column of graphs indicate travel time differences as a function of the takeoff-angles of the seismic rays.}
\label{Figure 1}
\end{center}
\end{figure}
\section{Least-squares method to find the best-fit solution}
\indent The predicted arrival time difference $ \Delta t^{i} $ at station i can be expressed as [Jordan \& Sverdrup, 1981]
\begin{align} \label{eq:1}
	\Delta t^{i} = \Delta t_{tra}^{i} + \Delta t_0, \ \  i = 1, 2, 3, \cdots, N,
	\end{align}
where N is the number of waveforms, and the predicted travel-time difference is defined by
\begin{align} \label{eq:1}
	\Delta t_{tra}^{i} = \mathbf{p}^{(i)} \cdot  \Delta  \mathbf{r},
	\end{align}
where $\mathbf{p}^{(i)} = (p_x^{i},p_y^{i},p_z^{i})$ is the i-th slowness vector, and $\Delta \mathbf{r} = (\Delta x, \Delta y, \Delta z)$ is the relative location between the two subevents. Correspondingly, 
\begin{align} \label{eq:1}
	p_x^{i} &= \frac{-\Delta x \sin(i_{0})_{i}\cos(az)_{i}}{v_p}, \\
	p_y^{i} &= \frac{-\Delta y \sin(i_{0})_{i}\cos(az)_{i}}{v_p}, \\
	p_z^{i} &= \frac{\Delta z \cos(i_{0})_{i}}{v_p}. 
	\end{align}
In this case, the forward problem can be expressed as
\begin{align} \label{eq:1}
	\mathbf{d} = \mathbf{G}\mathbf{m}. 
	\end{align}
Correspondingly, 
\begin{align} \label{eq:3}
	\mathbf{d} = 
	\begin{bmatrix}
	\Delta t^{(1)} \\
	\Delta t^{(2)} \\
	\Delta t^{(3)} \\
	\vdots \\
	\Delta t^{(N)} \\
	\end{bmatrix},
	\end{align}
	\begin{align} \label{eq:4}
	\mathbf{G} = 
	\begin{bmatrix}
	 p_x^{(1)} & p_y^{(1)}  & p_z^{(1)} & 1 \\
	 p_x^{(2)} & p_y^{(2)}  & p_z^{(2)} & 1 \\
	 p_x^{(3)} & p_y^{(3)}  & p_z^{(3)} & 1 \\
	\vdots & \vdots & \vdots & \vdots \\ 
	 p_x^{(N)} & p_y^{(N)}  & p_z^{(N)} & 1 \\
	\end{bmatrix}, 
	\end{align}
	and
	\begin{align} \label{eq:5}
	     \mathbf{m} = 
		\begin{bmatrix}
	 \Delta x \\
	 \Delta y\\
	 \Delta z\\
	\Delta t_0\\
	\end{bmatrix},
	\end{align}
	where $\Delta t^{(i)}$  can be estimated using waveform cross-correlation. \\
We assume the errors are Gaussian-distributed and solve equation 3.10 in a weighted least-squares sense: 
\begin{align} \label{eq:1}
	\mathbf{m}^{est} = (\mathbf{G}^{T}\mathbf{W}\mathbf{G})^{-1}\mathbf{G}^{T}\mathbf{W}\mathbf{d},
	\end{align}
 where $\mathbf{W}$ is the weighting matrix which is set up using a biweight function [Moesteller and Tukey, 1977; Waldhauser and Ellsworth, 2000]: 
 \begin{align} \label{eq:1}
	W_{i} = \text{max}^{2}\Bigg( 0,1-\Bigg( \frac{e_{i}}{\alpha\frac{\mathbf{e}_{MAD}}{\sigma_{MAD}}}\Bigg)^2 \Bigg).
\end{align}
$e_{i}$ represents the time residuals ($e_{i} = d_{i} - G_{ij}m^{est}_{j}$), $\mathbf{e}_{MAD} = \text{med}(\|e_i - \text{med}(\mathbf{e})\|)$ is the median absolute deviation from the median (MAD), $\sigma_{MAD} = 0.67449$ is the MAD for Gaussian noise, and $\alpha$ is a factor that defines the rejection level at $\alpha$ standard deviation. We set $\alpha$   to be 3. We set the maximum number of iterations to be 10. For the first run, the weighting matrix is set to be an $\text{N}\times \text{N}$ identity matrix. In the subsequent iterations,  the weighting matrix is determined by by time residuals distribution as in eq. 3.15.\\
 \section{Relative location error estimates}
 Even though we assume the errors are simply Gaussian-distributed, the relationship between the model parameters and the data can be complicated due to timing errors and errors due to an imprecise velocity model for the source region. It can cause the principle of error propagation to be inaccurate. Therefore, we apply a bootstrap method [Efron, 1982; Waldhauser \& Ellsworth, 2000] to estimate the location uncertainties empirically. For each run, we randomly choose a station from the original dataset with replacement until the resampled data size is equal to the original one. This process is repeated 10,000 times to ensure that the model parameters ($\Delta x, \Delta y, \Delta z, \Delta t_0$) are normally distributed. 
  \clearemptydoublepage

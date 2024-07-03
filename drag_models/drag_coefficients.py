def gamma1(x):
    return (1/(2*sqrt(pi)))*(exp(-x**2) + sqrt(pi)*x*(1+erf(x)))
def gamma2(x):
    return (1/(2*sqrt(pi)))*(x*exp(-x**2) + (sqrt(pi)/2)*(1 + 2*(x**2))*(1+erf(x)))
def C1(x):
    return (2*x**2 + 3)*iv(0,0.5*x**2) + (2*x**2 + 1)*iv(1,0.5*x**2)
def C2(x):
    return (x**2 + 1)*iv(0,0.5*x**2) + (x**2)*iv(1,0.5*x**2)
def drag_CLL_plane(v,R,temp,angle,sigma_t,sigma_n,T_s):
    s = v/sqrt(2*R*temp)
    v_w = sqrt(pi*R*T_s/2)
    g = cos(angle)
    c1 = sigma_t*(1 - g**2)*gamma1(g*s)
    c2 = (2-sigma_n)*g*gamma2(g*s)/s
    c3 = sigma_n*(v_w/v)*gamma1(g*s)*g
    return (2/s)*(c1+c2+c3)
def drag_CLL_sphere(v,R,temp,sigma_t,sigma_n,T_s):
    v_a = sqrt(2*R*temp);
    v_w = sqrt(pi*R*T_s/2);
    s = v/v_a
    p0 = (2 - sigma_n + sigma_t)/(2*(s**3))
    p1= (((4*(s**4)) + (4*(s**2)) - 1)*erf(s)/(2*s)) + ((2*s**2 + 1)*exp(-s**2)/sqrt(pi))
    p2 = (4/3)*sigma_n*(v_w/v) 
    return (p0*p1) + p2
def drag_CLL_cyl(v,R,temp,angle,sigma_t,sigma_n,T_s):
    ## angle = angle between cylindrical axis and velocity vector
    v_a = sqrt(2*R*temp);
    v_w = sqrt(pi*R*T_s/2);
    s = v/v_a
    g = sin(angle)
    mu = g*s
    p0 = (pi/2)*sigma_n*(v_w/v)*(g**2)
    p1 = (sqrt(pi)/(6*s))*(2*sigma_n - sigma_t - 4)*(g**2)*exp(-mu**2/2)*C1(mu)
    p2 = (sqrt(pi)*sigma_t/s)*(1 - g**2)*exp(-mu**2/2)*C2(mu)
    return p0 - p1 + p2
def drag_GRACE_CLL(v,temp,pitch,yaw,m_atm,nO,T_s,a,n_a,m_s): ## applicable for a flat-plate satellite model
    R_const = 8.31; R = (R_const/m_atm)*1000;
    pO = nO*k_B*temp; KL_CLL = 2.89*10**6;
    th = KL_CLL*pO/(1 + KL_CLL*pO);
    v_unit = (cos(-pitch)*cos(yaw),cos(-pitch)*sin(yaw),sin(-pitch))
    CD_CO = 0; CD_CL = 0;
    mu = m_atm/np.array(m_s) ## considers m_s is a list
    for i,n in enumerate(n_a):
        angle = arccos(dot(v_unit,n))
        alpha_n = (2*3.0*mu[i]/((1+mu[i])**2)) - 1
        if alpha_n < 0:
            alpha_n = 0
        sigma_n = 1 - sqrt(1-alpha_n)
        CD_CO = CD_CO + drag_CLL_plane(v,R,temp,angle,1,1,T_s)*a[i]
        CD_CL = CD_CL + drag_CLL_plane(v,R,temp,angle,1,sigma_n,T_s)*a[i]
    area = sat_area(pitch,yaw,a,n_a);
    C_D = (th*CD_CO + (1-th)*CD_CL)/area
    return C_D
def drag_DRIA_plane(v,R,temp,angle,alpha,T_s):
    s = v/sqrt(2*R*temp)
    g=cos(angle)
    mu = g*s
    P = exp(-mu**2)/s
    Q = 1 + 1/(2*(s**2))
    Z = 1 + erf(g*s)
    ratio = sqrt((2/3)*(1 + alpha*((3*R*T_s/(v**2))-1)))
    return (P/sqrt(pi)) + (g*Q*Z) + ((g/2)*ratio*((g*sqrt(pi)*Z) + P))
def drag_DRIA_sphere(v,R,temp,alpha,T_s):
    s = v/sqrt(2*R*temp)
    ratio = sqrt((2/3)*(1 + alpha*((3*R*T_s/(v**2))-1)))
    p0 = (2*s**2 + 1)*exp(-s**2)/(sqrt(pi) * s**3)
    p1 = (4*s**4 + 4*s**2 - 1)*erf(s)/(2 * s**4)
    p2 = (2*sqrt(pi)/3)*ratio
    return p0 + p1 + p2
def drag_DRIA_cyl(v,R,temp,angle,alpha,T_s):
    ## angle = angle between cylindrical axis and velocity vector
    s = v/sqrt(2*R*temp)
    g=sin(angle)
    mu = g*s
    ratio = sqrt((2/3)*(1 + alpha*((3*R*T_s/(v**2))-1)))
    p0 = s*sqrt(pi)*(g**2)*(1 + 1/(2 * s**2))*exp(-mu**2/2)*(iv(0,mu**2/2) + iv(1,mu**2/2))
    p1 = (sqrt(pi)/s)*exp(-mu**2/2)*iv(0,mu**2/2)
    p2 = (pi**1.5)*(g**2)*ratio/4
    return p0+p1+p2
def drag_GRACE_DRIA(v,temp,pitch,yaw,m_atm,nO,T_s,a,n_a,m_s):
    R_const = 8.31; R = (R_const/m_atm)*1000; 
    pO = nO*k_B*temp; KL_DRIA=1.44*10**6;
    th = KL_DRIA*pO/(1 + KL_DRIA*pO)
    v_unit = (cos(-pitch)*cos(yaw),cos(-pitch)*sin(yaw),sin(-pitch))
    mu = m_atm/np.array(m_s) ## considers m_s is a list
    CD_CO = 0; CD_CL = 0;
    for i,n in enumerate(n_a):
        angle = arccos(dot(v_unit,n))
        alpha = 3.0*mu[i]/((1+mu[i])**2)
        CD_CO = CD_CO + drag_DRIA_plane(v,R,temp,angle,1,T_s)*a[i]
        CD_CL = CD_CL + drag_DRIA_plane(v,R,temp,angle,alpha,T_s)*a[i]
    area = sat_area(pitch,yaw,a,n_a);
    C_D = (th*CD_CO + (1-th)*CD_CL)/area
    return C_D
def GRACE_area(pitch,yaw,a,n_a):
    area=0;
    v_unit = (cos(-pitch)*cos(yaw),cos(-pitch)*sin(yaw),sin(-pitch));
    for i,n in enumerate(n_a):
        angle = arccos(dot(v_unit,n))
        if angle < pi/2:
            area = area + a[i]*cos(angle)
    return area
def find_orbits(XYZ_TOD,start,stop):
    Z = XYZ_TOD['Z'][(XYZ_TOD.index >= start) & (XYZ_TOD.index <= stop)]
    idxz = np.where(np.diff(np.sign(Z)) > 0)[0]
    orbit_end = Z.index[idxz]
    time_avg = [orbit_end[i] + (orbit_end[i+1] - orbit_end[i])/2 for i in range(len(orbit_end)-1)]
    return orbit_end, time_avg
def itrf_to_gcrf(itrf_coords, times):

    """
    Convert ITRF coordinates to GCRF coordinates.

    Parameters:
    itrf_coords (numpy.ndarray): Array of ITRF coordinates (x, y, z) in meters.
    times (numpy.ndarray): Array of times in ISO format.

    Returns:
    numpy.ndarray: Array of GCRF coordinates (x, y, z) in meters.
    """
    # Initialize an empty list to hold the GCRF coordinates
    gcrf_coords = []

    # Loop over each set of ITRF coordinates and corresponding time
    for coords, time in zip(itrf_coords, times):
        # Convert ITRF coordinates to astropy CartesianRepresentation
        itrf_representation = CartesianRepresentation(*coords, unit=u.m)

        # Create ITRS frame
        itrs_frame = ITRS(itrf_representation, obstime=Time(time))

        # Convert to GCRS frame
        gcrs_frame = itrs_frame.transform_to(GCRS(obstime=Time(time)))

        # Extract GCRF coordinates and append to the list
        gcrf_coords.append(gcrs_frame.cartesian.xyz.to(u.m).value)

    # Convert the list of GCRF coordinates to a numpy array and return
    return np.array(gcrf_coords)
def get_dates(filename,start_year,stop_year):
    # Initialize empty lists to store the start dates and stop dates
    start_dates = []
    stop_dates = []

    # Open the file
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            # If the line contains a start date, extract the date
            if line.startswith('start'):
                # Extract the date from the line
                date = dt.datetime.strptime(line[19:-2], '(%Y,%m,%d)')
                # Add the date to the start dates list
                if date.year >= start_year and date.year <= stop_year:
                    start_dates.append(date)
            # If the line contains a stop date, extract the date
            elif line.startswith('stop'):
                # Extract the date from the line
                date = dt.datetime.strptime(line[18:-2], '(%Y,%m,%d)')
                # Add the date to the stop dates list
                if date.year >= start_year and date.year <= stop_year:
                    stop_dates.append(date)
    if len(start_dates) > len(stop_dates):
        start_dates = start_dates[:len(stop_dates)]

    # Return the start dates and stop dates
    return start_dates, stop_dates
def quaternion_to_euler(w, x, y, z):
    """
    Convert a quaternion into pitch, roll, and yaw angles.
    
    Parameters:
    w, x, y, z : float
        Quaternion components.
    
    Returns:
    pitch, roll, yaw : float
        Euler angles in radians.
    """
    
    # Roll (x-axis rotation)
    sinr_cosp = 2 * (w * x + y * z)
    cosr_cosp = 1 - 2 * (x * x + y * y)
    # atan2 ensures the angle is in the correct quadrant
    roll = math.atan2(sinr_cosp, cosr_cosp)
    
    # Pitch (y-axis rotation)
    sinp = 2 * (w * y - z * x)
    if abs(sinp) >= 1:
        # Handle out of range value for asin
        pitch = math.copysign(math.pi / 2, sinp)
    else:
        pitch = math.asin(sinp)
    
    # Yaw (z-axis rotation)
    siny_cosp = 2 * (w * z + x * y)
    cosy_cosp = 1 - 2 * (y * y + z * z)
    # atan2 ensures the angle is in the correct quadrant
    yaw = math.atan2(siny_cosp, cosy_cosp)
    
    return np.array([pitch, roll, yaw])
def CNOFS_area(v,a,n_a):
    area=0; r_probe = 0.06; r_boom = 0.0143;
    v_abs = np.sqrt(np.sum(v**2))
    v_unit = v/v_abs
    for i,n in enumerate(n_a):
        angle = arccos(dot(v_unit,n))
        if angle < pi/2:
            area = area + a[i]*cos(angle)
    pitch = -arcsin(v_unit[2]/v_abs)
    yaw = arctan(v_unit[1]/v_unit[0])
    area += 2*r_boom*(9.8-0.61)*cos((pi/4) - pitch)*2
    area += 2*r_boom*(9.8-0.61)*cos((pi/4) + pitch)*2
    area += 2*r_boom*18.02*cos(yaw)
    area += 6*(pi*(r_probe**2))
    return area
def drag_CNOFS_CLL(v,temp,m_atm,nO,T_s,a,n_a,m_s):
    r_probe = 0.06; r_boom = 0.0143;
    R_const = 8.31; R = (R_const/m_atm)*1000;
    pO = nO*k_B*temp; KL_CLL = 2.89*10**6;
    th = KL_CLL*pO/(1 + KL_CLL*pO);
    v_abs = np.sqrt(np.sum(v**2))
    v_unit = v/v_abs
    pitch = -arcsin(v_unit[2]/v_abs)
    yaw = arctan(v_unit[1]/v_unit[0])
    CD_CO = 0; CD_CL = 0;
    mu = m_atm/m_s;
    alpha_n = (2*3.0*mu/((1+mu)**2)) - 1
    if alpha_n < 0:
        alpha_n = 0
    sigma_n = 1 - sqrt(1-alpha_n)
    for i,n in enumerate(n_a):
        angle = arccos(dot(v_unit,n))
        CD_CO = CD_CO + drag_CLL_plane(v_abs,R,temp,angle,1,1,T_s)*a[i]
        CD_CL = CD_CL + drag_CLL_plane(v_abs,R,temp,angle,1,sigma_n,T_s)*a[i]

    CD_CO += 2*drag_CLL_cyl(v_abs,R,temp,pi/4 + pitch,1,1,T_s)*2*r_boom*(9.8-0.61)
    CD_CO += 2*drag_CLL_cyl(v_abs,R,temp,pi/4 - pitch,1,1,T_s)*2*r_boom*(9.8-0.61)

    CD_CL += 2*drag_CLL_cyl(v_abs,R,temp,pi/4 + pitch,1,sigma_n,T_s)*2*r_boom*(9.8-0.61)
    CD_CL += 2*drag_CLL_cyl(v_abs,R,temp,pi/4 - pitch,1,sigma_n,T_s)*2*r_boom*(9.8-0.61)
    
    CD_CO += drag_CLL_cyl(v_abs,R,temp,pi/2 - yaw,1,1,T_s)*2*r_boom*18.02
    CD_CL += drag_CLL_cyl(v_abs,R,temp,pi/2 - yaw,1,sigma_n,T_s)*2*r_boom*18.02
    
    CD_CO += (6*drag_CLL_sphere(v_abs,R,temp,1,1,T_s)*(pi*(r_probe**2)))
    CD_CL += (6*drag_CLL_sphere(v_abs,R,temp,1,sigma_n,T_s)*(pi*(r_probe**2)))
    area = CNOFS_area(v,a,n_a);
    C_D = (th*CD_CO + (1-th)*CD_CL)/area
    return C_D
def drag_CNOFS_DRIA(v,temp,m_atm,nO,T_s,a,n_a,m_s):
    r_probe = 0.06; r_boom = 0.0143;
    R_const = 8.31; R = (R_const/m_atm)*1000;
    pO = nO*k_B*temp; KL_DRIA = 1.44*10**6;
    th = KL_DRIA*pO/(1 + KL_DRIA*pO);
    v_abs = np.sqrt(np.sum(v**2))
    v_unit = v/v_abs
    pitch = -arcsin(v_unit[2]/v_abs)
    yaw = arctan(v_unit[1]/v_unit[0])
    CD_CO = 0; CD_CL = 0;
    mu = m_atm/m_s;
    alpha_s = 3.0*mu/((1+mu)**2)
    
    for i,n in enumerate(n_a):
        angle = arccos(dot(v_unit,n))
        CD_CO = CD_CO + drag_DRIA_plane(v_abs,R,temp,angle,1,T_s)*a[i]
        CD_CL = CD_CL + drag_DRIA_plane(v_abs,R,temp,angle,alpha_s,T_s)*a[i]

    CD_CO += 2*drag_DRIA_cyl(v_abs,R,temp,pi/4 + pitch,1,T_s)*2*r_boom*(9.8-0.61)
    CD_CO += 2*drag_DRIA_cyl(v_abs,R,temp,pi/4 - pitch,1,T_s)*2*r_boom*(9.8-0.61)

    CD_CL += 2*drag_DRIA_cyl(v_abs,R,temp,pi/4 + pitch,alpha_s,T_s)*2*r_boom*(9.8-0.61)
    CD_CL += 2*drag_DRIA_cyl(v_abs,R,temp,pi/4 - pitch,alpha_s,T_s)*2*r_boom*(9.8-0.61)
    
    CD_CO += drag_DRIA_cyl(v_abs,R,temp,pi/2 - yaw,1,T_s)*2*r_boom*18.02
    CD_CL += drag_DRIA_cyl(v_abs,R,temp,pi/2 - yaw,alpha_s,T_s)*2*r_boom*18.02
    
    CD_CO += (6*drag_DRIA_sphere(v_abs,R,temp,1,T_s)*(pi*(r_probe**2)))
    CD_CL += (6*drag_DRIA_sphere(v_abs,R,temp,alpha_s,T_s)*(pi*(r_probe**2)))
    area = CNOFS_area(v,a,n_a);
    C_D = (th*CD_CO + (1-th)*CD_CL)/area
    return C_D




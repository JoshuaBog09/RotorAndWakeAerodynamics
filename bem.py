

def inflow_velocity(a:float, fv_free:float):
    vp = (1-a)*v_free
    return vp

def lift(c:float, rho:float, cl:float, vp:float):
    lift = 0.5*rho*vp**2*c*cl
    return lift

def drag(c:float, rho:float, cd:float, vp:float):
    drag = 0.5*rho*vp**2*c*cd
    return lift

def force_azimut(lift:float, drag:float, twist:float):
    f_azi = lift*np.sin(twist)-drag*np.cos(twist)
    return f_azi

def force_axi(lift:float, drag:float, twist:float):
    f_azi = lift*np.cos(twist)+drag*np.sin(twist)
    return f_azi
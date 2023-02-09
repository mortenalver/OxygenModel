package fishmodel.sim;

/**
* Created by malv on 12.06.15.
*/ // Variable simulation settings. Class used as struct:
public class SimulationInputs implements Cloneable {
    public double sinkingSpeed;
    public int centerX, centerY;
    public double feedingRateMult = 0;
    public double x_vel = 0.05,
            y_vel = 0;
    public double x_windSpeed=0, y_windspeed=0;
    public int airSpeed = 30;
    public int spreaderType = 0;
    public boolean tiltUp = true;

    @Override
    public Object clone() {
        SimulationInputs si = new SimulationInputs();
        si.sinkingSpeed = sinkingSpeed;
        si.centerX = centerX;
        si.centerY = centerY;
        si.feedingRateMult = feedingRateMult;
        si.x_vel = x_vel;
        si.y_vel = y_vel;
        si.x_windSpeed = x_windSpeed;
        si.y_windspeed = y_windspeed;
        si.airSpeed = airSpeed;
        si.spreaderType = spreaderType;
        si.tiltUp = tiltUp;
        return si;
    }
}

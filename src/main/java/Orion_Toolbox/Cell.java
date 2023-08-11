
package Orion_Toolbox;

import java.util.HashMap;
import mcib3d.image3d.ImageHandler;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.measurements.MeasureIntensity;


/**
 * @author orion-cirb
 */
public class Cell {
    
    private Object3DInt cell;
    private Object3DInt nucleus;
    private Object3DInt cytoplasm;
    public HashMap<String, Double> params;
    
    public Cell(Object3DInt cell, Object3DInt nucleus, Object3DInt cytoplasm) {
        this.cell = cell;
        this.nucleus = nucleus;
        this.cytoplasm = cytoplasm;
        this.params = new HashMap<>();
    }
    
    public void fillVolumes(double pixelVol) {
        params.put("cellVol", cell.size() * pixelVol);
        params.put("nucleusVol", nucleus.size() * pixelVol);
        params.put("cytoplasmVol", cytoplasm.size() * pixelVol);
    }
    
    public void fillIntensities(ImageHandler imh) {
        params.put("cellInt", new MeasureIntensity(cell, imh).getValueMeasurement(MeasureIntensity.INTENSITY_AVG));
        params.put("nucleusInt", new MeasureIntensity(nucleus, imh).getValueMeasurement(MeasureIntensity.INTENSITY_AVG));
        params.put("cytoplasmInt", new MeasureIntensity(cytoplasm, imh).getValueMeasurement(MeasureIntensity.INTENSITY_AVG));
    }
}

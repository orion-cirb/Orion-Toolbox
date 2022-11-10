package Tools;

import Stardist.StarDist2D;
import Cellpose.CellposeTaskSettings;
import Cellpose.CellposeSegmentImgPlusAdvanced;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.ZProjector;
import ij.process.ImageProcessor;
import ij.util.ThreadUtil;
import java.awt.Font;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom2.BoundingBox;
import mcib3d.geom2.Object3DComputation;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Object3DPlane;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.VoxelInt;
import mcib3d.geom2.measurements.Measure2Colocalisation;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.geom2.measurementsPopulation.PairObjects3DInt;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.spatial.analysis.SpatialStatistics;
import mcib3d.spatial.descriptors.F_Function;
import mcib3d.spatial.descriptors.SpatialDescriptor;
import mcib3d.spatial.sampler.SpatialModel;
import mcib3d.spatial.sampler.SpatialRandomHardCore;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;


/**
 * @author Orion-CIRB
 */
public class Tools {

    private final File stardistModelsPath = new File(IJ.getDirectory("imagej")+File.separator+"models");
    public Calibration cal;
    
     // Omnipose
    private final String omniposeEnvDirPath = (IJ.isLinux()) ? "/opt/miniconda3/envs/omnipose" : System.getProperty("user.home")+"\\miniconda3\\envs\\OmniPose";
    
    // Cellpose
    private final String cellposeEnvDirPath = (IJ.isLinux()) ? "/opt/miniconda3/envs/cellpose" : System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose";
    
    
    /**
     * Check that needed modules are installed
     * @param className
     * @param PluginName
     * @return 
     */
    public boolean checkInstalledModules(String className, String PluginName) {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass(className);
        } catch (ClassNotFoundException e) {
            IJ.showMessage("Error", PluginName+" not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    /**
     * Flush and close an image
     * @param img
     */
    public void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Find images extension
     * @param imagesFolder
     * @return 
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        File[] files = imagesFolder.listFiles();
        for (File file: files) {
            if(file.isFile()) {
                String fileExt = FilenameUtils.getExtension(file.getName());
                switch (fileExt) {
                   case "nd" :
                       ext = fileExt;
                       break;
                   case "nd2" :
                       ext = fileExt;
                       break;
                    case "czi" :
                       ext = fileExt;
                       break;
                    case "lif"  :
                        ext = fileExt;
                        break;
                    case "ics2" :
                        ext = fileExt;
                        break;
                    case "tif" :
                        ext = fileExt;
                        break;
                    case "tiff" :
                        ext = fileExt;
                        break;
                }
            } else if (file.isDirectory() && !file.getName().equals("Results")) {
                ext = findImageType(file);
                if (! ext.equals(""))
                    break;
            }
        }
        return(ext);
    }
     
    /**
     * Check that required StarDist models are present in Fiji models folder
     * @param stardistModel
     * @return 
     */
    public boolean checkStarDistModels(String stardistModel) {
        FilenameFilter filter = (dir, name) -> name.endsWith(".zip");
        File[] modelList = stardistModelsPath.listFiles(filter);
        int index = ArrayUtils.indexOf(modelList, new File(stardistModelsPath+File.separator+stardistModel));
        if (index == -1) {
            IJ.showMessage("Error", stardistModel + " StarDist model not found, please add it in Fiji models folder");
            return false;
        }
        return true;
    }
    
    /**
     * Find images in folder
     * @param imagesFolder
     * @param imageExt
     * @return 
     */
    public ArrayList<String> findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in "+imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
    /**
     * Find image calibration
     * @param meta
     * @return 
     */
    public Calibration findImageCalib(IMetadata meta) {
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
        return(cal);
    }
    
    
     /**
     * Find channels name
     * @param imageName
     * @param meta
     * @param reader
     * @return 
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                }
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                }
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelFluor(0, n);
                break;
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;    
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        return(channels);         
    }
    
    /**
     * Do Z projection
     * @param img
     * @param param
     * @return 
     */
    public ImagePlus doZProjection(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
    
     /**
     * Find background image intensity:
     * Z projection over min intensity
     * return Mean + std of intensity
     * @param img
     * @return 
     */
    public double findMeanStdBackground(ImagePlus img) {
      double[] bg = new double[2];
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      bg[0] = imp.getStatistics().mean;
      bg[1] = imp.getStatistics().stdDev;
      System.out.println("Background (mean +- std of the min projection) = " + bg[0] + " +- " + bg[1]);
      flush_close(imgProj);
      return(bg[0]+bg[1]);
    }
    
    /**
     * Find background image intensity:
     * Z projection over min intensity
     * return median of intensity
     * @param img
     * @return 
     */
    public double findMedianBackground(ImagePlus img) {
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      double bg = imp.getStatistics().mean;
      System.out.println("Background (median of the min projection) = " + bg);
      flush_close(imgProj);
      return(bg);
    }
    
     /**
     * Clear out side roi
     * @param img
     * @param roi
     */
    public void clearOutSide(ImagePlus img, Roi roi) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), Roi.FREEROI);
        poly.setLocation(0, 0);
        for (int n = 1; n <= img.getNSlices(); n++) {
            ImageProcessor ip = img.getImageStack().getProcessor(n);
            ip.setRoi(poly);
            ip.setBackgroundValue(0);
            ip.setColor(0);
            ip.fillOutside(poly);
        }
        img.updateAndDraw();
    }
   
    
     /**
     * Label object
     * @param obj
     * @param img 
     * @param fontSize 
     */
    public void labelObject(Object3DInt obj, ImagePlus img, int fontSize) {
        if (IJ.isMacOSX())
            fontSize *= 3;
        
        BoundingBox bb = obj.getBoundingBox();
        int z = bb.zmin + 1;
        int x = bb.xmin;
        int y = bb.ymin;
        img.setSlice(z);
        ImageProcessor ip = img.getProcessor();
        ip.setFont(new Font("SansSerif", Font.PLAIN, fontSize));
        ip.setColor(255);
        ip.drawString(Integer.toString((int)obj.getLabel()), x, y);
        img.updateAndDraw();
    }
    
    /**
     * compute local thickness
     * @param img
     * @param inverse
     * @return 
    **/
    public ImageFloat localThickness3D (ImagePlus img, boolean inverse) {
        IJ.showStatus("Computing distance map...");
        img.setCalibration(cal);
        ImageFloat edt = new EDT().run(ImageHandler.wrap(img), 0,inverse, ThreadUtil.getNbCpus());
        return(edt);
    }
    
    /**
     * Remove object with only one plan
     * @param pop
     */
    public void popFilterOneZ(Objects3DIntPopulation pop) {
        pop.getObjects3DInt().removeIf(p -> (p.getObject3DPlanes().size() == 1));
        pop.resetLabels();
    }
    
    /**
     * Remove object with size < min and size > max
     * @param pop
     * @param min
     * @param max
     */
    public void popFilterSize(Objects3DIntPopulation pop, double min, double max) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
        pop.resetLabels();
    }
    
    /**
     * Take only object above an intensity Threshold
     * @param imgLabels
     * @param imgRaw
     * @param intensityThresh
     * @return 
     */
    public ImagePlus filterDetectionsByIntensity(ImagePlus imgLabels, ImagePlus imgRaw, double intensityThresh) {
        double background = findMeanStdBackground(imgRaw);
        Objects3DIntPopulation popIn = new Objects3DIntPopulation(ImageHandler.wrap(imgLabels.getStack()));
        Objects3DIntPopulation popOut = new Objects3DIntPopulation();
        ImageHandler imhRaw = ImageHandler.wrap(imgRaw);
        for (Object3DInt obj: popIn.getObjects3DInt()) {
            if (new MeasureIntensity(obj, imhRaw).getValueMeasurement(MeasureIntensity.INTENSITY_AVG) > intensityThresh*background) {
                popOut.addObject(obj);
            }
        }
        return(popOut.drawImage().getImagePlus());
    }
    
     /**
     * Find volume of objects  
     * @param dotsPop
     * @return vol
     */
    
    public double findPopVolume (Objects3DIntPopulation dotsPop) {
        IJ.showStatus("Findind object's volume");
        List<Double[]> results = dotsPop.getMeasurementsList(new MeasureVolume().getNamesMeasurement());
        double sum = results.stream().map(arr -> arr[1]).reduce(0.0, Double::sum);
        return(sum);
    }
    
     /**
     * Find intensity of objects  
     * @param dotsPop
     * @param img
     * @return intensity
     */
    
    public double findPopIntensity (Objects3DIntPopulation dotsPop, ImagePlus img) {
        IJ.showStatus("Findind object's intensity");
        ImageHandler imh = ImageHandler.wrap(img);
        double sumInt = 0;
        for(Object3DInt obj : dotsPop.getObjects3DInt()) {
            MeasureIntensity intMes = new MeasureIntensity(obj, imh);
            sumInt +=  intMes.getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
        }
        return(sumInt);
    }
    
    /**
     * Find objects in pop colocalized with obj
     * %coloc of object in pop > pourc
     * @param obj
     * @param pop
     * @param pourc
     * @return 
     */
    public Objects3DIntPopulation findColocObjects(Object3DInt obj, Objects3DIntPopulation pop, double pourc) {
        Objects3DIntPopulation colocPop = new Objects3DIntPopulation();
        if (pop.getNbObjects() > 0) {
            for (Object3DInt objPop: pop.getObjects3DInt()) {
                double colocPour = new Measure2Colocalisation(objPop, obj).getValue(Measure2Colocalisation.COLOC_PC);
                if (colocPour >= pourc) {
                    colocPop.addObject(objPop);
                }
            }
        }
        colocPop.resetLabels();
        return(colocPop);
    }
    
     /**
     * Find coloc between pop1 and pop2
     * set label of colocalized object in Id object
     * @param pop1
     * @param pop2
     * @param pourc
     * @return number of pop objects colocalized with pop1
     * @throws java.io.IOException
     */
    public int findNumberColocPop (Objects3DIntPopulation pop1, Objects3DIntPopulation pop2, double pourc) throws IOException {
        if (pop1.getNbObjects() == 0 && pop2.getNbObjects() == 0) 
            return(0);
        MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(pop1, pop2);
        AtomicInteger ai = new AtomicInteger(0);
        pop1.getObjects3DInt().forEach(obj1 -> {
            float obj1label = obj1.getLabel();
            List<PairObjects3DInt> list = coloc.getPairsObject1(obj1label, true);
            if (!list.isEmpty()) {
                list.forEach(p -> {
                    Object3DInt obj2 = p.getObject3D2();
                    if (p.getPairValue() > obj2.size()*pourc) {
                        obj2.setIdObject(obj1label);
                        ai.incrementAndGet();
                    }
                });
            }
        });
        return(ai.get());
    } 
    
    
    /**
     * Find coloc objects in pop1 colocalized with pop2
     * @param pop1
     * @param pop2
     * @param pourc
     * @return 
     */
    public Objects3DIntPopulation findColocPop (Objects3DIntPopulation pop1, Objects3DIntPopulation pop2, double pourc) {
        Objects3DIntPopulation colocPop = new Objects3DIntPopulation();
        if (pop1.getNbObjects() > 0 && pop2.getNbObjects() > 0) {
            MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(pop1, pop2);
            pop1.getObjects3DInt().forEach(obj1 -> {
                List<PairObjects3DInt> list = coloc.getPairsObject1(obj1.getLabel(), true);
                if (!list.isEmpty()) {
                    list.forEach(p -> {
                        Object3DInt obj2 = p.getObject3D2();
                        if (p.getPairValue() > obj2.size()*pourc) {
                            colocPop.addObject(obj2);
                        }
                    });
                }
            });
        }
        return(colocPop);
    }
    
    /**
     * Return dilated object restriced to image borders
     * @param img
     * @param obj
     * @param dilSize
     * @return 
     */
    public Object3DInt dilateObj(ImagePlus img, Object3DInt obj, double dilSize) {
        Object3DInt objDil = new Object3DComputation(obj).getObjectDilated((float)(dilSize/cal.pixelWidth), (float)(dilSize/cal.pixelHeight), 
                (float)(dilSize/cal.pixelDepth));
        // check if object go outside image
        BoundingBox bbox = objDil.getBoundingBox();
        BoundingBox imgBbox = new BoundingBox(ImageHandler.wrap(img));
        int[] box = {imgBbox.xmin, imgBbox.xmax, imgBbox.ymin, imgBbox.ymax, imgBbox.zmin, imgBbox.zmax};
        if (bbox.xmin < 0 || bbox.xmax > imgBbox.xmax || bbox.ymin < 0 || bbox.ymax > imgBbox.ymax
                || bbox.zmin < 0 || bbox.zmax > imgBbox.zmax) {
            Object3DInt objDilImg = new Object3DInt();
            for (Object3DPlane p : objDil.getObject3DPlanes()) {
                for (VoxelInt v : p.getVoxels()) {
                    if (v.isInsideBoundingBox(box))
                        objDilImg.addVoxel(v);
                }
            }
            return(objDilImg);
        }
        else
            return(objDil);
    }
    
    /**
     * Apply StarDist 2D slice by slice
     * Label detections in 3D
     * @param img
     * @param factor
     * @param resize
     * @param blockRad
     * @param stardistModel
     * @param stardistProbThresh
     * @param stardistOverlayThresh
     * @return objects population
     * @throws java.io.IOException
     */
    public Objects3DIntPopulation stardistObjectsPop(ImagePlus img, float factor, boolean resize, int blockRad, String stardistModel,
            double stardistProbThresh, double stardistOverlayThresh) throws IOException {
        Object syncObject = new Object();
        double stardistPercentileBottom = 0.2;
        double stardistPercentileTop = 99.8;
        String stardistOutput = "Label Image";
        
        // Resize image to be in a StarDist-friendly scale
        int imgWidth = img.getWidth();
        int imgHeight = img.getHeight();
        ImagePlus imgIn = null;
        if (resize) {
            if (factor > 1) {
                imgIn = img.resize((int)(imgWidth*factor), (int)(imgHeight*factor), 1, "bicubic");
            }
            else 
                imgIn = new ImagePlus("", img.getProcessor().resize((int)(imgWidth*factor), (int)(imgHeight*factor), false));
            resize = true;
        }
        else
            imgIn = new Duplicator().run(img);
        
        if (blockRad != 0)
             // Remove outliers
             IJ.run(imgIn, "Remove Outliers", "block_radius_x="+blockRad+" block_radius_y="+blockRad+" standard_deviations=1 stack");

        // StarDist
        File starDistModelFile = new File(stardistModelsPath+File.separator+stardistModel);
        StarDist2D star = new StarDist2D(syncObject, starDistModelFile);
        star.loadInput(imgIn);
        star.setParams(stardistPercentileBottom, stardistPercentileTop, stardistProbThresh, stardistOverlayThresh, stardistOutput);
        star.run();
        flush_close(imgIn);

        // Label detections in 3D
        ImagePlus imgOut = (resize) ? star.getLabelImagePlus().resize(imgWidth, imgHeight, 1, "none") : star.getLabelImagePlus();       
        ImagePlus imgLabels = star.associateLabels();
        imgLabels.setCalibration(cal); 
        flush_close(imgOut);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgLabels));      
        flush_close(imgLabels);
       return(pop);
    }
    
     /**
    * Detect Cells with CellPose
     * @param img
     * @param cellposeModel
     * @param cellposeDiameter
     * @param cellposeMaskThreshold
     * @param cellposeFlowThreshold
     * @param factor
     * @param resize
     * @param useGPU
     * @return 
    */
    public Objects3DIntPopulation cellposeDetection(ImagePlus img, String cellposeModel, int cellposeDiameter, 
            double cellposeMaskThreshold, double cellposeFlowThreshold, int factor, boolean resize, boolean useGPU){
        // Resize to be in a cellpose-friendly scale
        ImagePlus imgIn = null;
        int imgWidth = img.getWidth();
        int imgHeight = img.getHeight();
        if (resize) {
            if (factor > 1) {
                imgIn = img.resize((int)(imgWidth*factor), (int)(imgHeight*factor), 1, "bicubic");
            }
            else 
                imgIn = new ImagePlus("", img.getProcessor().resize((int)(imgWidth*factor), (int)(imgHeight*factor), false));
            resize = true;
        }
        else
            imgIn = new Duplicator().run(img);
        
        imgIn.setCalibration(cal);
        
        // Set Cellpose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(cellposeModel, 1, cellposeDiameter, cellposeEnvDirPath);
        settings.setCellProbTh(cellposeMaskThreshold);
        settings.setFlowTh(cellposeFlowThreshold);
        settings.useGpu(useGPU);
        
        // Run Omnipose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgIn);
        ImagePlus imgOut = cellpose.run();
        
        ImageProcessor imgOutProc = (resize) ? imgOut.getProcessor().resize(imgWidth, imgHeight, false) : imgOut.getProcessor();
        imgOut = new ImagePlus("", imgOutProc);
        imgOut.setCalibration(cal);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgOut));
        
        // Close images
        flush_close(imgIn);
        flush_close(imgOut);
        return(pop);
    }
    
    
    /**
    * Detect bacteria with Omnipose
     * @param img
     * @param omniposeModel
     * @param omniposeDiameter
     * @param omniposeMaskThreshold
     * @param omniposeFlowThreshold
     * @param factor
     * @param resize
     * @param useGPU
     * @return 
    */
    public Objects3DIntPopulation omniposeDetection(ImagePlus img, String omniposeModel, int omniposeDiameter, 
            double omniposeMaskThreshold, double omniposeFlowThreshold, int factor, boolean resize, boolean useGPU){
        // Resize to be in a Omnipose-friendly scale
        ImagePlus imgIn = null;
        int imgWidth = img.getWidth();
        int imgHeight = img.getHeight();
        if (resize) {
            if (factor > 1) {
                imgIn = img.resize((int)(imgWidth*factor), (int)(imgHeight*factor), 1, "bicubic");
            }
            else 
                imgIn = new ImagePlus("", img.getProcessor().resize((int)(imgWidth*factor), (int)(imgHeight*factor), false));
            resize = true;
        }
        else
            imgIn = new Duplicator().run(img);
        
        imgIn.setCalibration(cal);
        
        // Set Omnipose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(omniposeModel, 1, omniposeDiameter, omniposeEnvDirPath);
        settings.setVersion("0.7");
        settings.setCluster(true);
        settings.setOmni(true);
        settings.useMxNet(false);
        settings.setCellProbTh(omniposeMaskThreshold);
        settings.setFlowTh(omniposeFlowThreshold);
        settings.useGpu(useGPU);
        
        // Run Omnipose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgIn);
        ImagePlus imgOut = cellpose.run();
        
        ImageProcessor imgOutProc = (resize) ? imgOut.getProcessor().resize(imgWidth, imgHeight, false) : imgOut.getProcessor();
        imgOut = new ImagePlus("", imgOutProc);
        imgOut.setCalibration(cal);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgOut));
        
        // Close images
        flush_close(imgIn);
        flush_close(imgOut);
        return(pop);
    }
    
    /**
     * Compute F-function-related Spatial Distribution Index of foci population in a nucleus
     * https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000853
     */ 
    public Double computeSdiF(Objects3DIntPopulation fociInt, Object3DInt nucInt, ImagePlus img) {
        // Convert Object3DInt & Objects3DIntPopulation objects into Object3D & Objects3DPopulation objects
        ImageHandler imhNuc = ImageHandler.wrap(img).createSameDimensions();
        nucInt.drawObject(imhNuc, 1);
        Object3D nuc = new Objects3DPopulation(imhNuc).getObject(0);
        ImageHandler imhFoci = ImageHandler.wrap(img).createSameDimensions();
        fociInt.drawInImage(imhFoci);
        Objects3DPopulation foci = new Objects3DPopulation(imhFoci);
        // Define spatial descriptor and model
        SpatialDescriptor spatialDesc = new F_Function(2500, nuc); // nb of points used to compute the F-function        
        SpatialModel spatialModel = new SpatialRandomHardCore(foci.getNbObjects(), 0.8/cal.pixelWidth, nuc); // average diameter of a spot in pixels
        SpatialStatistics spatialStatistics = new SpatialStatistics(spatialDesc, spatialModel, 100, foci); // nb of samples (randomized organizations simulated to compare with the spatial organization of the spots)
        spatialStatistics.setEnvelope(0.05); // 2.5-97.5% envelope error
        spatialStatistics.setVerbose(false);
        return(spatialStatistics.getSdi());
    }
    

    
    
}

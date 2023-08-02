package Orion_Tools;

import Orion_Tools.StardistOrion.StarDist2D;
import Orion_Tools.Cellpose.CellposeTaskSettings;
import Orion_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import features.TubenessProcessor;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.ZProjector;
import ij.plugin.filter.BackgroundSubtracter;
import ij.process.ImageProcessor;
import java.awt.Font;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.DoubleAccumulator;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom2.BoundingBox;
import mcib3d.geom2.Object3DComputation;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Object3DIntLabelImage;
import mcib3d.geom2.Object3DPlane;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.VoxelInt;
import mcib3d.geom2.measurements.Measure2Colocalisation;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.geom2.measurementsPopulation.PairObjects3DInt;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.processing.BinaryMorpho;
import mcib3d.spatial.descriptors.G_Function;
import mcib3d.spatial.descriptors.F_Function;
import mcib3d.spatial.descriptors.SpatialDescriptor;
import mcib3d.spatial.sampler.SpatialModel;
import mcib3d.spatial.sampler.SpatialRandomHardCore;
import mcib3d.utils.ThreadUtil;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import net.haesleinhuepf.clijx.imagej2.ImageJ2Tubeness;
import net.haesleinhuepf.clijx.bonej.BoneJSkeletonize3D;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;


/**
 * @author Orion-CIRB
 */
public class Tools {
    
    private final File stardistModelsPath = new File(IJ.getDirectory("imagej")+File.separator+"models");
    private Object syncObject = new Object();
        
     // Omnipose
    private final String omniposeEnvDirPath = (IJ.isLinux()) ? "/opt/miniconda3/envs/omnipose" : System.getProperty("user.home")+"\\miniconda3\\envs\\OmniPose";
    
    // Cellpose
    public final String cellposeEnvDir = IJ.isWindows()? System.getProperty("user.home")+File.separator+"miniconda3"+File.separator+"envs"+File.separator+"CellPose" : "/opt/miniconda3/envs/cellpose";
    private final String cellposeModelsPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\.cellpose\\models\\" : System.getProperty("user.home")+"/.cellpose/models/"; 
    
    public Calibration cal;
    private final CLIJ2 clij2 = CLIJ2.getInstance();

    
    /**
     * Display a message in the ImageJ console and status bar
     * @param log
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
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
     * Find stardist models in Fiji models folder
     * 
     * @return 
    */
    public String[] findStardistModels() {
        FilenameFilter filter = (dir, name) -> name.endsWith(".zip");
        File[] modelList = stardistModelsPath.listFiles(filter);
        String[] modelFiles = new String[modelList.length];
        for (int i = 0; i < modelList.length; i++)
            modelFiles[i] = modelList[i].getName();
        return modelFiles;
    }

    
    
    /**
    ********************************************************
     * Clij filters
    ********************************************************
     */
    
    
    /**
     * Test if GPU
     * @return 
     */
    public boolean isGPU() {
        String gpuName = clij2.getGPUName();
        return(gpuName.isEmpty());
    }
    
    /**
     * 3D Median filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    public ImagePlus median3D_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img); 
       ClearCLBuffer imgCLMed = clij2.create(imgCL);
       clij2.median3DBox(imgCL, imgCLMed, sizeXY, sizeXY, sizeZ);
       clij2.release(imgCL);
       ImagePlus imgMed = clij2.pull(imgCLMed);
        clij2.release(imgCLMed);
       return(imgMed);
    } 
    
    /**
     * 2D Median filter using CLIJ2
     * @param img
     * @param sizeXY
     * @return 
     */ 
    public ImagePlus median2D_filter(ImagePlus img, double sizeXY) {
       ClearCLBuffer imgCL = clij2.push(img); 
       ClearCLBuffer imgCLMed = clij2.create(imgCL);
       clij2.median2DBox(imgCL, imgCLMed, sizeXY, sizeXY);
       clij2.release(imgCL);
       ImagePlus imgMed = clij2.pull(imgCLMed);
        clij2.release(imgCLMed);
       return(imgMed);
    } 

    /**
     * 3D Gaussian filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    public ImagePlus gaussian3D_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img); 
       ClearCLBuffer imgCLGauss = clij2.create(imgCL);
       clij2.gaussianBlur3D(imgCL, imgCLGauss, sizeXY, sizeXY, sizeZ);
       clij2.release(imgCL);
       ImagePlus imgGauss = clij2.pull(imgCLGauss);
       clij2.release(imgCLGauss);
       return(imgGauss);
    } 
    
    
     /**
     * Difference of Gaussians 
     * Using CLIJ2
     * @param img
     * @param size1
     * @param size2
     * @return imgGauss
     */ 
    public ImagePlus DOG(ImagePlus img, double size1, double size2) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian3D(imgCL, imgCLDOG, size1, size1, size1, size2, size2, size2);
        clij2.release(imgCL);
        ImagePlus imgDOG = clij2.pull(imgCLDOG);
        clij2.release(imgCLDOG);
        return(imgDOG);
    }
    
    /**
     * Threshold 
     * USING CLIJ2
     * @param img
     * @param thMed
     * @return 
     */
    public ImagePlus threshold(ImagePlus img, String thMed) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        clij2.release(imgCL);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCLBin);
        return(imgBin);
    }
    
     /**
     * Laplace of Gaussian filter using CLIJ2
     * @param img
     * @param sizeXYZ
     * @return 
     */ 
    public ImagePlus LOG_filter(ImagePlus img, double sizeXYZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLDOG = clij2.create(imgCL);
       clij2.gaussianBlur3D(imgCL, imgCLDOG, sizeXYZ, sizeXYZ, sizeXYZ);
       clij2.release(imgCL);
       ClearCLBuffer imgCLLOG = clij2.create(imgCL);
       clij2.laplaceSphere(imgCLDOG, imgCLLOG);
       clij2.release(imgCLDOG);
       ImagePlus imgLOG = clij2.pull(imgCLLOG);
       clij2.release(imgCLLOG);
       return(imgLOG);
    }  
    
    
    /**
     * 3D Max filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    public ImagePlus max3D_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMax = clij2.create(imgCL);
       clij2.maximum3DBox(imgCL, imgCLMax, sizeXY, sizeXY, sizeZ);
       clij2.release(imgCL);
       return(clij2.pull(imgCLMax));
    } 
    
    /**
     * 3D Min filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    public ImagePlus min3D_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMin = clij2.create(imgCL);
       clij2.minimum3DBox(imgCL, imgCLMin, sizeXY, sizeXY, sizeZ);
       clij2.release(imgCL);
       return(clij2.pull(imgCLMin));
    } 
    
    /**
     * 3D Opening filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    public ImagePlus open3D_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMin = clij2.create(imgCL);
       clij2.minimum3DBox(imgCL, imgCLMin, sizeXY, sizeXY, sizeZ);
       clij2.maximum3DBox(imgCLMin, imgCL, sizeXY, sizeXY, sizeZ);
       clij2.release(imgCLMin);
       return(clij2.pull(imgCL));
    }
    
    /**
     * 3D closing filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    public ImagePlus close3D_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMax = clij2.create(imgCL);
       clij2.maximum3DBox(imgCLMax, imgCL, sizeXY, sizeXY, sizeZ);
       clij2.minimum3DBox(imgCL, imgCLMax, sizeXY, sizeXY, sizeZ);
       clij2.release(imgCLMax);
       return(clij2.pull(imgCL));
    }

    /**
     * tubeness plane by plane multithreads
     */
    public ImagePlus tubeness(ImagePlus img, float sigma) {
        ImagePlus[] stack = new ImagePlus[img.getNSlices()];
        final AtomicInteger ai = new AtomicInteger(0);
        final int n_cpus = ThreadUtil.getNbCpus();
        Thread[] threads = ThreadUtil.createThreadArray(n_cpus);
        int nbZ = img.getNSlices();
        final int dec = (int) Math.ceil((double) nbZ / (double) n_cpus);
        
        for (int ithread = 0; ithread < threads.length; ithread++) {
            threads[ithread] = new Thread() {
                @Override
                public void run() {
                    for (int k = ai.getAndIncrement(); k < n_cpus; k = ai.getAndIncrement()) {
                        for (int z = dec * k; ((z < (dec * (k + 1))) && (z < nbZ)); z++) {
                            ImagePlus slice = new Duplicator().run(img, z+1, z+1);
                            stack[z] = new TubenessProcessor(sigma, true).generateImage(slice);
                        }
                    }
                }
            };
        }
        ThreadUtil.startAndJoin(threads);
        ImageStack imgS = new ImageStack(img.getWidth(), img.getHeight());
        for (int i = 0; i < nbZ; i++)
           imgS.addSlice(stack[i].getProcessor());
        return(new ImagePlus("", imgS));
    }

    /**
     * subtract background plane by plane multithreads
     */
    public ImagePlus subtractBackground(ImagePlus img, float radius) {
        final AtomicInteger ai = new AtomicInteger(0);
        final int n_cpus = ThreadUtil.getNbCpus();
        Thread[] threads = ThreadUtil.createThreadArray(n_cpus);
        int nbZ = img.getNSlices();
        final int dec = (int) Math.ceil((double) nbZ / (double) n_cpus);
        ImagePlus imgBG = img.duplicate();
        for (int ithread = 0; ithread < threads.length; ithread++) {
            threads[ithread] = new Thread() {
                @Override
                public void run() {
                    for (int k = ai.getAndIncrement(); k < n_cpus; k = ai.getAndIncrement()) {
                        for (int z = dec * k; ((z < (dec * (k + 1))) && (z < nbZ)); z++) {
                            imgBG.setZ(z+1);
                            ImageProcessor ip = imgBG.getProcessor();
                            new BackgroundSubtracter().rollingBallBackground​(ip, radius, false, false, false, false, false);
                        }
                    }
                }
            };
        }
        ThreadUtil.startAndJoin(threads);
        return(imgBG);
    }
    
    /**
     * Clij2 skeletonize 3D
     */
    public ImagePlus skeletonize3D(ImagePlus img) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLSkel = clij2.create(imgCL);
        BoneJSkeletonize3D skel = new BoneJSkeletonize3D();
        skel.bonejSkeletonize3D(clij2,imgCL, imgCLSkel);
        clij2.release(imgCL);
        ImagePlus imgSkel = clij2.pull(imgCLSkel);
        clij2.release(imgCLSkel);
        return(imgSkel);
    }

   
      /**
     * Morphological operators with clij
     * @param obj
     * @param op
     * @param rad (pixels)
     * @return return null object if touch border
     */
    
    public Object3DInt getMorphologicalObject3D(Object3DInt obj, ImagePlus img, int op, int rad) {
        int ext = (op == BinaryMorpho.MORPHO_DILATE) ? rad + 1 : 1;
        ImageHandler labelImage = new Object3DIntLabelImage(obj).getCroppedLabelImage(ext, ext, 0, 1, false);
        ImagePlus imgCrop = labelImage.getImagePlus();
        ImagePlus imgSeg = null;
        switch (op) {
            case BinaryMorpho.MORPHO_DILATE :
                imgSeg = max3D_filter(imgCrop, rad, 0);
                break;
            case BinaryMorpho.MORPHO_ERODE :
                imgSeg = min3D_filter(imgCrop, rad, 0);
                break;
        }
        ImageHandler segImage2 = ImageHandler.wrap(imgSeg);
        segImage2.setOffset(labelImage);
        segImage2.setCalibration(cal);
        Object3DInt objMorpho = new Object3DInt(segImage2);
        if ((op == BinaryMorpho.MORPHO_DILATE) && (new Object3DComputation(objMorpho).touchBorders(ImageHandler.wrap(img), false)))
            objMorpho = null;
        else
            objMorpho.setLabel(obj.getLabel());
        return objMorpho;
    }
    
    
    /**
    ********************************************************
     * File tools
    ********************************************************
     */
    
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
                    case "ics" :
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
            } else if (file.isDirectory() && !file.getName().contains("Results")) {
                ext = findImageType(file);
                if (! ext.equals(""))
                    break;
            }
        }
        return(ext);
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
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).getValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).getValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
        return(cal);
    }
    
    /**
     * Find channels name and None to end of list
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
        String[] channels = new String[chs+1];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).equals("")) ? Integer.toString(n) : meta.getChannelName(0, n);
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).equals("")) ? Integer.toString(n) : meta.getChannelName(0, n);
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
                    channels[n] = (meta.getChannelFluor(0, n).equals("")) ? Integer.toString(n) : meta.getChannelFluor(0, n);
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).toString();
                break;    
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).toString();
                break; 
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        channels[chs] = "None";
        return(channels);     
    }
    
    
    /**
     ********************************************************
     * Processing tools
     ******************************************************** 
     */
    
    
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
     * Find background image intensity
     * with/without roi
     * Z projection over min intensity + read stats (median/mean) intensity
     * @param img
     * @param roi
     * @param method
     * @return 
     */
    public double findBackground(ImagePlus img, Roi roi, String method) {
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      if (roi != null) {
          roi.setLocation(0, 0);
          imp.setRoi(roi);
      }
      double bg = (method.equals("median")) ? imp.getStatistics().median : imp.getStatistics().mean;
      System.out.println("Background = " + bg);
      flush_close(imgProj);
      return(bg);
    }
    
    /**
     * Find background image intensity
     * with/without roi
     * Z projection over min intensity + read stats (median/mean) intensity
     * @param img
     */
    public double findStdDevBackground(ImagePlus img, Roi roi) {
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      if (roi != null) {
          roi.setLocation(0, 0);
          imp.setRoi(roi);
      }
      double bg = imp.getStatistics().stdDev;
      System.out.println("Background = " + bg);
      flush_close(imgProj);
      return(bg);
    }
    
    
    /**
     * Auto find background from scroolling roi
     * @param img
     * @param roiBgSize
     * @param method
     * @return 
     */
    public RoiBg findRoiBackgroundAuto(ImagePlus img, int roiBgSize, String method) {
        // scroll image and measure bg intensity in roi 
        // take roi lower intensity
        
        ArrayList<RoiBg> intBgFound = new ArrayList<>();
        for (int x = 0; x < img.getWidth() - roiBgSize; x += roiBgSize) {
            for (int y = 0; y < img.getHeight() - roiBgSize; y += roiBgSize) {
                Roi roi = new Roi(x, y, roiBgSize, roiBgSize);
                double bg = findBackground(img, roi, method);
                intBgFound.add(new RoiBg(roi, bg));
            }
        }
        img.deleteRoi();
        // sort RoiBg on bg value
        intBgFound.sort(Comparator.comparing(RoiBg::getBgInt));
        // Find lower value
        RoiBg roiBg = intBgFound.get(0);
        
        int roiCenterX = (int)(roiBg.getRoi().getBounds().x+(roiBgSize/2));
        int roiCenterY = (int)(roiBg.getRoi().getBounds().y+(roiBgSize/2));
        System.out.println("Roi auto background("+method+") found = "+roiBg.getBgInt()+" center x = "+roiCenterX+", y = "+roiCenterY);
        return(roiBg);
    }
    
    
     /**
     * Clear out side roi
     * @param img
     * @param roi
     */
    public void clearOutSide(ImagePlus img, Roi roi) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), roi.getType());
        roi.setLocation(0, 0);
        for (int n = 1; n <= img.getNSlices(); n++) {
            ImageProcessor ip = img.getImageStack().getProcessor(n);
            ip.setRoi(poly);
            ip.setBackgroundValue(0);
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
    public void labelObjectCenter(Object3DInt obj, ImagePlus img, int fontSize) {
        if (IJ.isMacOSX())
            fontSize *= 3;
        
        BoundingBox bb = obj.getBoundingBox();
        int z = (int)(bb.zmin + 0.5*(bb.zmax - bb.zmin));
        int x = (int)(bb.xmin + 0.5*(bb.xmax - bb.xmin));
        int y = (int)(bb.ymin + 0.5*(bb.ymax - bb.ymin));
        img.setSlice(z); //z+1
        ImageProcessor ip = img.getProcessor();
        ip.setFont(new Font("SansSerif", Font.PLAIN, fontSize));
        ip.setColor(255);
        ip.drawString(Integer.toString((int)obj.getLabel()), x, y);
        img.updateAndDraw();
    }
    
     /**
     * Label object
     * @param obj
     * @param img 
     * @param fontSize 
     */
    public void labelObjectLeftTop(Object3DInt obj, ImagePlus img, int fontSize) {
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
     * Compute distance map or inverse distance map
     * with clij2
     * @param img
     * @param inverse
     * @return 
     */
    public ImagePlus localThickness3D(ImagePlus img, boolean inverse) {
        IJ.showStatus("Computing distance map...");
        ImagePlus imgDup = new Duplicator().run(img);
        IJ.run(imgDup, "8-bit", "");
        if (inverse)
            IJ.run(imgDup, "Invert", "stack");
        ClearCLBuffer imgCL = clij2.push(imgDup);
        flush_close(imgDup);
        ClearCLBuffer imgCLMap = clij2.create(imgCL);
        clij2.distanceMap(imgCL, imgCLMap);
        clij2.release(imgCL);
        ImagePlus imgMap = clij2.pull(imgCLMap);
        clij2.release(imgCLMap);
        imgMap.setCalibration(cal);
        return(imgMap);
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
     * Remove object touching border image
     */
    public void removeTouchingBorder(Objects3DIntPopulation pop, ImagePlus img) {
        ImageHandler imh = ImageHandler.wrap(img);
        pop.getObjects3DInt().removeIf(p -> (new Object3DComputation(p).touchBorders(imh, false)));
        pop.resetLabels();
    }
    
    
      /**
     * Remove object with size < min and size > max
     * @param pop
     * @param min
     * @param max
     */
    public void popFilterSize(Objects3DIntPopulation pop, double min, double max) {
        pop.setVoxelSizeXY(cal.pixelWidth);
        pop.setVoxelSizeZ(cal.pixelDepth);
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
        pop.resetLabels();
    }
    
    /**
     * Remove objects in population with intensity < intTh
     * @param pop
     * @param img
     * @param intTh 
     */
    public void intensityFilter(Objects3DIntPopulation pop, ImagePlus img, double intTh) {
        ImageHandler imh = ImageHandler.wrap(img);
        pop.getObjects3DInt().removeIf(p -> (new MeasureIntensity(p, imh).getValueMeasurement(MeasureIntensity.INTENSITY_MAX) < intTh));
        pop.resetLabels();
    }
    
    /**
     * Find sum vessel volume
     */
    private double findPopVolume(Objects3DIntPopulation pop) {
        DoubleAccumulator sumVol = new DoubleAccumulator(Double::sum,0.d);
        pop.getObjects3DInt().parallelStream().forEach(obj -> { 
            sumVol.accumulate(new MeasureVolume(obj).getVolumeUnit());
        });
        return(sumVol.doubleValue());
    }
    
     /**
     * Find sum intensity of objects  
     * @param dotsPop
     * @param img
     * @return intensity
     */
    private double findPopIntensity(Objects3DIntPopulation pop, ImagePlus img) {
        ImageHandler imh = ImageHandler.wrap(img);
        DoubleAccumulator sumInt = new DoubleAccumulator(Double::sum,0.d);
        pop.getObjects3DInt().parallelStream().forEach(obj -> { 
            sumInt.accumulate(new MeasureIntensity(obj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM));
        });
        return(sumInt.doubleValue());
    }
    
    /**
     * Find objects in pop colocalized with obj
     * %coloc of object in pop object > pourc
     * @param obj
     * @param pop
     * @param pourc
     * @return 
     */
    public Objects3DIntPopulation findColocObjects(Object3DInt obj, Objects3DIntPopulation pop, double pourc) {
        Objects3DIntPopulation colocPop = new Objects3DIntPopulation();
        if (pop.getNbObjects() > 0) {
            pop.getObjects3DInt().parallelStream().filter(obj2 -> {
                double colocPour = new Measure2Colocalisation(obj, obj2).getValue(Measure2Colocalisation.COLOC_PC);
                return(colocPour >= pourc);
            }).forEach(colocPop::addObject);
        colocPop.resetLabels();
        colocPop.setVoxelSizeXY(cal.pixelWidth);
        colocPop.setVoxelSizeZ(cal.pixelDepth);
    }
        return(colocPop);
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
            pop1.getObjects3DInt().stream().forEach(obj1 -> {
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
        colocPop.resetLabels();
        return(colocPop);
    }

   
    
    /** Dilate objects
     * Return croped object to image borders or remove object touching border 
     * A revoir !!!!!
     * @param img
     * @param pop
     * @param dilSize
     * @return 
     */
    public Object3DInt restrictedObjToImg(ImagePlus img, Objects3DIntPopulation pop, int dilSize, boolean removeObj) {
        // check if object go outside image
        Object3DInt objDil = new Object3DInt();
        for (Object3DInt obj : pop.getObjects3DInt()){
            BoundingBox bbox = obj.getBoundingBox();
            BoundingBox imgBbox = new BoundingBox(ImageHandler.wrap(img));
            int[] box = {imgBbox.xmin, imgBbox.xmax, imgBbox.ymin, imgBbox.ymax, imgBbox.zmin, imgBbox.zmax};
            if (bbox.xmin < 0 || bbox.xmax > imgBbox.xmax || bbox.ymin < 0 || bbox.ymax > imgBbox.ymax
                    || bbox.zmin < 0 || bbox.zmax > imgBbox.zmax || !removeObj) {
                for (Object3DPlane p : obj.getObject3DPlanes()) {
                    for (VoxelInt v : p.getVoxels()) {
                        if (v.isInsideBoundingBox(box))
                            objDil.addVoxel(v);
                    }
                }
            }
            else {
                removeTouchingBorder(pop, img);
            }
        }
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
        String method = (factor > 1) ? "bicubic" : "none";
        ImagePlus imgIn = (resize) ? img.resize((int)(imgWidth*factor), (int)(imgHeight*factor), 1, method) : new Duplicator().run(img);
        
        if (blockRad != 0)
             // Remove outliers
            IJ.run(imgIn, "Remove Outliers...", "radius="+blockRad+" threshold=1 which=Bright stack");

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

        // Resize image to be in a StarDist-friendly scale
        int imgWidth = img.getWidth();
        int imgHeight = img.getHeight();
        String method = (factor > 1) ? "bicubic" : "none";
        ImagePlus imgIn = (resize) ? img.resize((int)(imgWidth*factor), (int)(imgHeight*factor), 1, method) : new Duplicator().run(img);
        String model = cellposeModelsPath+cellposeModel;
        // Set Cellpose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(model, 1, cellposeDiameter, cellposeEnvDir);
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
        
        // Resize image to be in a StarDist-friendly scale
        int imgWidth = img.getWidth();
        int imgHeight = img.getHeight();
        String method = (factor > 1) ? "bicubic" : "none";
        ImagePlus imgIn = (resize) ? img.resize((int)(imgWidth*factor), (int)(imgHeight*factor), 1, method) : new Duplicator().run(img);
        
        // Set Omnipose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(omniposeModel, 1, omniposeDiameter, omniposeEnvDirPath);
        settings.setVersion("0.7");
        settings.setCluster(true);
        settings.setOmni(true);
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
     * Compute G-function-related Spatial Distribution Index of cells population in a Obj/ROI
     * @param popInt
     * @param roiInt
     * @param img
     * @param distHardCore
     * @param numRandomSamples
     * @param plotName
     * @return
     */
    public double[] computeSdiG(Objects3DIntPopulation popInt, Object3DInt roiInt, ImagePlus img, double distHardCore, int numRandomSamples, String plotName) {
        // Convert Object3DInt & Objects3DIntPopulation objects into Object3D & Objects3DPopulation objects
        ImageHandler imhRoi = ImageHandler.wrap(img).createSameDimensions();
        roiInt.drawObject(imhRoi, 1);
        Object3D roi = new Objects3DPopulation(imhRoi).getObject(0);
        ImageHandler imhPop = ImageHandler.wrap(img).createSameDimensions();
        popInt.drawInImage(imhPop);
        Objects3DPopulation pop = new Objects3DPopulation(imhPop);
        
        // Define spatial descriptor and model
        SpatialDescriptor spatialDesc = new G_Function();     
        SpatialModel spatialModel = new SpatialRandomHardCore(pop.getNbObjects(), distHardCore, roi); // average diameter of a cell in pixels
        SpatialStatistics spatialStatistics = new SpatialStatistics(spatialDesc, spatialModel, numRandomSamples, pop); // nb of samples (randomized organizations simulated to compare with the spatial organization of the cells)
        spatialStatistics.setEnvelope(0.05); // 2.5-97.5% envelope error
        spatialStatistics.setVerbose(false);
        double sdiG = spatialStatistics.getSdi();
        double area = spatialStatistics.getAreaCurve();
        
        Plot plotG = spatialStatistics.getPlot();
        plotG.draw();
        plotG.addLabel(0.1, 0.1, "SDI = " + String.format("%.3f", sdiG));
        ImagePlus imgPlot = plotG.getImagePlus();
        FileSaver plotSave = new FileSaver(imgPlot);
        plotSave.saveAsTiff(plotName);
        flush_close(imgPlot); 
        double[] results = {sdiG, area};
        return(results);
    }
 /**
     * Compute F-function-related Spatial Distribution Index of cells population in a Obj/ROI
     * @param popInt
     * @param roiInt
     * @param img
     * @param distHardCore
     * @param numRandomSamples
     * @param plotName
     * @return
     */
    public double[] computeSdiF(Objects3DIntPopulation popInt, Object3DInt roiInt, ImagePlus img, double distHardCore, int numRandomSamples, String plotName) {
        // Convert Object3DInt & Objects3DIntPopulation objects into Object3D & Objects3DPopulation objects
        ImageHandler imhRoi = ImageHandler.wrap(img).createSameDimensions();
        roiInt.drawObject(imhRoi, 1);
        Object3D roi = new Objects3DPopulation(imhRoi).getObject(0);
        ImageHandler imhPop = ImageHandler.wrap(img).createSameDimensions();
        popInt.drawInImage(imhPop);
        Objects3DPopulation pop = new Objects3DPopulation(imhPop);
        
        // Define spatial descriptor and model
        SpatialDescriptor spatialDesc = new F_Function(2500, roi); // nb of points used to compute the F-function    
        SpatialModel spatialModel = new SpatialRandomHardCore(pop.getNbObjects(), distHardCore, roi); // average diameter of a cell in pixels
        SpatialStatistics spatialStatistics = new SpatialStatistics(spatialDesc, spatialModel, numRandomSamples, pop); // nb of samples (randomized organizations simulated to compare with the spatial organization of the cells)
        spatialStatistics.setEnvelope(0.05); // 2.5-97.5% envelope error
        spatialStatistics.setVerbose(false);
        double sdiG = spatialStatistics.getSdi();
        double area = spatialStatistics.getAreaCurve();
        
        Plot plotG = spatialStatistics.getPlot();
        plotG.draw();
        plotG.addLabel(0.1, 0.1, "SDI = " + String.format("%.3f", sdiG));
        ImagePlus imgPlot = plotG.getImagePlus();
        FileSaver plotSave = new FileSaver(imgPlot);
        plotSave.saveAsTiff(plotName);
        flush_close(imgPlot); 
        double[] results = {sdiG, area};
        return(results);
    }
    
    
    
}

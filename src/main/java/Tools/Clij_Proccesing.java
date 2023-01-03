/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package Tools;

import ij.ImagePlus;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;

/**
 *
 * @author phm
 */
public class Clij_Proccesing {
    private final CLIJ2 clij2 = CLIJ2.getInstance();
    
    /**
     * Median filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    public ClearCLBuffer median_filter(ClearCLBuffer imgCL, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCLMed = clij2.create(imgCL);
       clij2.median3DBox(imgCL, imgCLMed, sizeXY, sizeXY, sizeZ);
       return(imgCLMed);
    }  
    
    /**
     * Test if GPU
     * @return 
     */
    public boolean isGPU() {
        String gpuName = clij2.getGPUName();
        return(gpuName.isEmpty());
    }
    
     /**
     * Difference of Gaussians 
     * Using CLIJ2
     * @param imgCL
     * @param size1
     * @param size2
     * @return imgGauss
     */ 
    public ClearCLBuffer DOG(ClearCLBuffer imgCL, double size1, double size2) {
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian3D(imgCL, imgCLDOG, size1, size1, size1, size2, size2, size2);
        return(imgCLDOG);
    }
    
    /**
     * Threshold 
     * USING CLIJ2
     * @param imgCL
     * @param thMed
     */
    public ClearCLBuffer threshold(ClearCLBuffer imgCL, String thMed) {
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        return(imgCLBin);
    }
    
     /**
     * Laplace of Gaussian filter using CLIJ2
     * @param img
     * @param sizeXYZ
     * @return 
     */ 
    public ImagePlus  LOG_filter(ImagePlus img, double sizeXYZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLDOG = clij2.create(imgCL);
       clij2.gaussianBlur3D(imgCL, imgCLDOG, sizeXYZ, sizeXYZ, sizeXYZ);
       clij2.release(imgCL);
       ClearCLBuffer imgCLLOG = clij2.create(imgCL);
       clij2.laplaceBox(imgCLDOG, imgCLLOG);
       clij2.release(imgCLDOG);
       ImagePlus imgLOG = clij2.pull(imgCLLOG);
       clij2.release(imgCLLOG);
       return(imgLOG);
    }  
    
    
}

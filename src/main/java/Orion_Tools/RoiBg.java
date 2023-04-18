/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Orion_Tools;

import ij.gui.Roi;

/**
 *
 * @author phm
 */
public class RoiBg {
    private Roi roi;
    private double bgInt;
    
	
	public RoiBg(Roi roi, double bgInt) {
            this.roi = roi;
            this.bgInt = bgInt;
	}
        
        public void setRoi(Roi roi) {
		this.roi = roi;
	}
        
        public void setBgInt(double bgInt) {
		this.bgInt = bgInt;
	}
        
        public Roi getRoi() {
            return roi;
        }
        
        public double getBgInt() {
            return bgInt;
        }

}


// Maximum Inscribed Circle plugin for OpenJUMP - Java version
//
// Copyright (C) 2020 Burak Beyhan, Cüneyt Güler and Hidayet Tağa
// 
// This plugin calculates Maximum Inscribed Circle (MIC) for polygon features. For further
// explanations please refer to the following papers if you use this plugin in your studies;
//
// Beyhan, B., Güler, C., Tağa, H. (2020) An algorithm for maximum inscribed circle based on Voronoi 
// diagrams and geometrical properties. Journal of Geographical Systems, 22, 391–418.
// https://doi.org/10.1007/s10109-020-00325-3
//
// This program is free software; you can redistribute it and/or modify it under the terms of 
// the GNU General Public License as published by the Free Software Foundation; either version 2 
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see http://www.gnu.org/licenses


package MIC.tools;
import java.awt.Color;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Iterator;
import javax.swing.JFrame;
import com.vividsolutions.jump.feature.AttributeType;
import com.vividsolutions.jump.feature.Feature;
import com.vividsolutions.jump.feature.FeatureCollection;
import com.vividsolutions.jump.feature.FeatureSchema;
import com.vividsolutions.jump.feature.FeatureDataset;
import com.vividsolutions.jump.feature.BasicFeature;
import com.vividsolutions.jump.workbench.model.Layer;
import com.vividsolutions.jump.workbench.model.LayerManager;
import com.vividsolutions.jump.workbench.model.LayerStyleUtil;
import com.vividsolutions.jump.workbench.plugin.PlugInContext;
import com.vividsolutions.jump.workbench.ui.GUIUtil;
import com.vividsolutions.jump.workbench.ui.renderer.style.BasicStyle;
import com.vividsolutions.jump.workbench.ui.task.TaskMonitorDialog;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.operation.buffer.BufferOp;


public class MICEngine {

  private Layer katmanim;
  private boolean usefield;
  private boolean addMIC;
  private boolean addCIM;
  private String attribim;
  private String cunitim;
  private String runitim;
  private int precdis;
  private Double csr;

  public MICEngine() {
  }
 
  public void aktarim1(Layer katman){katmanim = katman;}
  public void aktarim2(String st1){attribim = st1;}
  public void aktarim3(String st2){cunitim = st2;}
  public void aktarim4(boolean bo1){addMIC = bo1;}
  public void aktarim5(boolean bo2){addCIM = bo2;}
  public void aktarim6(String st3){runitim = st3;}
  public void aktarim7(boolean bo3){usefield = bo3;} 
  public void aktarim8(Double do1){csr = do1;}
  public void aktarim9(int in1){precdis = in1;} 

  public void execute(PlugInContext context) throws Exception {

    final LayerManager layerManager = context.getLayerManager();
    String activefilename = katmanim.getName();
    Layer actualLayer = katmanim;	
    FeatureCollection fc = actualLayer.getFeatureCollectionWrapper().getWrappee();
    FeatureSchema fs = fc.getFeatureSchema();
   
    String strob = new String(new char[precdis]).replace("\0", "#"); // creation of a decimal format for displaying the results in the attribute table
    DecimalFormat df = new DecimalFormat("#."+strob); // in line with the precision defined for this purpose
    DecimalFormatSymbols custom = new DecimalFormatSymbols();
    custom.setDecimalSeparator('.');
    df.setDecimalFormatSymbols(custom);
    
    double CBU; // unit of measurement that will be used in displaying the perimeter and radius of MIC in the attribute table
    if (cunitim == "millimeter (mm)") CBU = 0.001;
    else if (cunitim == "centimeter (cm)") CBU = 0.01;
    else if (cunitim == "decimeter (dm)") CBU = 0.1;
    else if (cunitim == "meter (m)") CBU = 1;
    else if (cunitim == "hectometer (hm)") CBU = 100;
    else CBU = 1000;
	double ABU; // unit of measurement that will be used in displaying the area of MIC in the attribute table
    if (runitim == "square mm (mm²)") ABU = 0.000001;
    else if (runitim == "square cm (cm²)")  ABU = 0.0001;
    else if (runitim == "square dm (dm²)")  ABU = 0.01;
    else if (runitim == "square m (m²)")  ABU = 1;
    else if (runitim == "decare (daa)")  ABU = 1000;
    else if (runitim == "hectare (ha)")  ABU = 10000;
    else ABU = 1000000;
    
   	int onlem = fs.getAttributeCount(); // number of existing fields in the attribute table

    FeatureSchema micSchema = new FeatureSchema(); // new feature schema for maximum inscribed circles
    micSchema.addAttribute(fs.getAttributeName(fs.getGeometryIndex()), AttributeType.GEOMETRY);   
    FeatureSchema picSchema = new FeatureSchema(); // new feature schema for the centers of maximum inscribed circles
    picSchema.addAttribute(fs.getAttributeName(fs.getGeometryIndex()), AttributeType.GEOMETRY);   
    
    Object typim = AttributeType.INTEGER; // if no field is selected for feature id to be used in enumeration (IDM) of the resulting features, it is created ... 
    if (usefield) typim = fs.getAttributeType(attribim);
    
    // adding parameters as new fields to the attribute tables
    if (!fs.hasAttribute("micperim")) fs.addAttribute("micperim", AttributeType.DOUBLE); // in the existing attribute table of the input layer
    if (!fs.hasAttribute("micarea")) fs.addAttribute("micarea", AttributeType.DOUBLE); // the parameters calculated for MIC may already exist
    if (!fs.hasAttribute("micradius")) fs.addAttribute("micradius", AttributeType.DOUBLE); // thus for each of these fields, this is checked
    if (!fs.hasAttribute("micX")) fs.addAttribute("micX", AttributeType.DOUBLE);
    if (!fs.hasAttribute("micY")) fs.addAttribute("micY", AttributeType.DOUBLE);
    micSchema.addAttribute("IDM", (AttributeType) typim); // parameters for MIC layer are added to attribute table
    micSchema.addAttribute("micperim", AttributeType.DOUBLE);
    micSchema.addAttribute("micarea", AttributeType.DOUBLE);
    micSchema.addAttribute("micradius", AttributeType.DOUBLE);
    micSchema.addAttribute("micX", AttributeType.DOUBLE);
    micSchema.addAttribute("micY", AttributeType.DOUBLE);
    picSchema.addAttribute("IDM", (AttributeType) typim); // parameters for MIC's center layer are added to attribute table
    picSchema.addAttribute("micradius", AttributeType.DOUBLE);
    picSchema.addAttribute("micX", AttributeType.DOUBLE);
    picSchema.addAttribute("micY", AttributeType.DOUBLE);
    
    FeatureCollection resultmic = new FeatureDataset(micSchema); // creation of feature collection for MICs
    FeatureCollection resultpic = new FeatureDataset(picSchema); // creation of feature collection for centers of MICs
    FeatureCollection newDataset = new FeatureDataset(fs); // creation of new feature collection for the input layer
 
    final JFrame desktop = (JFrame) context.getWorkbenchFrame();
    final TaskMonitorDialog progressDialog = new TaskMonitorDialog(desktop, null);
    progressDialog.setTitle("Processing Polygon Features");
    progressDialog.addComponentListener(new ComponentAdapter() {
      public void componentShown(ComponentEvent e) {
        new Thread(new Runnable() {
          public void run() {
            try {
            	
              // Beginning the iteration for the polygon features in the input layer  
              int sayar = 1; // initial value for IDM if no field is selected for feature id to be used in enumeration
              for (Iterator iter = fc.iterator(); iter.hasNext();) {
            	 Feature element = (Feature) iter.next();
            	 Object value = sayar;
            	 progressDialog.report(value+". Polygon is being processed ...");
            	 if (usefield) value = element.getAttribute(attribim);
            	 sayar++;
  	    
            	 BasicFeature nf = new BasicFeature(fs); // creation of a basic feature to store the changes that will be made 
            	 for (int i = 0 ; i < (onlem) ; i++) { // adding existing fields to the new basic feature
            		 nf.setAttribute(i, element.getAttribute(i));
            	 }
      
            	 final GeometryFactory gf = new GeometryFactory(); // creating a geometry factory for subsequent calculations
            	 Feature feature = new BasicFeature(micSchema);
            	 Feature peature = new BasicFeature(picSchema);
            	 Polygon geom = (Polygon) element.getGeometry(); // for processing polygon with holes ... 
      
            	 double[] gelo = MICFunction.getMIC(geom, csr, 2500);
            	 double micX = gelo[0];
            	 double micY = gelo[1];
            	 double radius = gelo[2];
            	 double micPeri = 2*Math.PI*radius;
            	 double micArea = Math.PI*radius*radius;
            	 com.vividsolutions.jts.geom.Point micnt = gf.createPoint(new Coordinate(micX, micY)); // creating centroid of MIC
            	 Geometry ebid = BufferOp.bufferOp(micnt, radius, 128); // creating MIC
	  
            	 // Writing MIC parameters to the attribute tables of the input layer
            	 nf.setGeometry(geom);
            	 nf.setAttribute("micperim", Double.parseDouble(df.format(micPeri/CBU)));
            	 nf.setAttribute("micarea", Double.parseDouble(df.format(micArea/ABU)));
            	 nf.setAttribute("micradius", Double.parseDouble(df.format(radius/CBU)));
            	 nf.setAttribute("micX", Double.parseDouble(df.format(micX)));
            	 nf.setAttribute("micY", Double.parseDouble(df.format(micY)));
            	 
                 if (addMIC) { //  writing MIC parameters to the attribute table of the layer created for Maximum Inscribed Circles
                	feature.setGeometry(ebid);
                	feature.setAttribute("IDM", value);
                	feature.setAttribute("micperim", Double.parseDouble(df.format(micPeri/CBU)));
                	feature.setAttribute("micarea", Double.parseDouble(df.format(micArea/ABU)));
                	feature.setAttribute("micradius", Double.parseDouble(df.format(radius/CBU)));
                	feature.setAttribute("micX", Double.parseDouble(df.format(micX)));
                	feature.setAttribute("micY", Double.parseDouble(df.format(micY)));
               	 	resultmic.add(feature);
                 }
                 
                 if (addCIM) { // writing MIC parameters to the attribute table of the layer created for the Centers of Maximum Inscribed Circles
                	peature.setGeometry(micnt);
                	peature.setAttribute("IDM", value);
                	peature.setAttribute("micradius", Double.parseDouble(df.format(radius/CBU)));
                	peature.setAttribute("micX", Double.parseDouble(df.format(micX)));
                	peature.setAttribute("micY", Double.parseDouble(df.format(micY)));
               	 	resultpic.add(peature);
                 }
	  
            	 newDataset.add(nf);
            	 
              } // end of the loop created for iteration of the polygon features in the input layer 
	
              
              actualLayer.setFeatureCollection(newDataset); // updating the input layer
              if (addMIC || addCIM) LayerStyleUtil.setLinearStyle(actualLayer, Color.black, 1, 0); // setting the style of the input layer
              actualLayer.fireAppearanceChanged();
              
              if (addMIC) { // if it is confirmed, the layer created for Maximum Inscribed Circles is added to the map canvas with a predefined style.
          	     Layer MIClyr = layerManager.addLayer("MIC Analysis for "+activefilename, "Maximum Inscribed Circles", resultmic);
          	     BasicStyle BS =  MIClyr.getBasicStyle();
          	     BS.setFillColor(Color.cyan);
          	     BS.setLineColor(Color.blue);
          	     BS.setAlpha(164);
          	     MIClyr.addStyle(BS);
          	     MIClyr.fireAppearanceChanged();
              }
              if (addCIM) { // if it is confirmed, the layer created for Centers of Maximum Inscribed Circles is added to the map canvas with a predefined style.
            	  Layer CIMlyr = layerManager.addLayer("MIC Analysis for "+activefilename, "Centers of Maximum Inscribed Circles", resultpic);
            	  BasicStyle CV =  CIMlyr.getBasicStyle();
            	  CV.setLineColor(Color.red);
            	  CV.setFillColor(Color.red);
            	  CV.setLineWidth(1);
            	  CIMlyr.addStyle(CV);
            	  CIMlyr.fireAppearanceChanged();
              }
              
            } catch (Exception e) {
            } finally {
                   progressDialog.setVisible(false);
                   progressDialog.dispose();
            }
          }
        }).start();
      }
    });
	
    GUIUtil.centreOnWindow(progressDialog);
    progressDialog.setVisible(true);
   
    return;
  }
}
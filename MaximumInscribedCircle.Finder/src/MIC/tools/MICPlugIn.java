
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
import java.util.ArrayList;
import java.util.List;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.DefaultComboBoxModel;
import javax.swing.ImageIcon;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import com.vividsolutions.jump.feature.AttributeType;
import com.vividsolutions.jump.feature.FeatureSchema;
import com.vividsolutions.jump.workbench.WorkbenchContext;
import com.vividsolutions.jump.workbench.model.Layer;
import com.vividsolutions.jump.workbench.plugin.AbstractPlugIn;
import com.vividsolutions.jump.workbench.plugin.EnableCheckFactory;
import com.vividsolutions.jump.workbench.plugin.MultiEnableCheck;
import com.vividsolutions.jump.workbench.plugin.PlugInContext;
import com.vividsolutions.jump.workbench.ui.GUIUtil;
import com.vividsolutions.jump.workbench.ui.MultiInputDialog;
import com.vividsolutions.jump.workbench.ui.images.IconLoader;
import com.vividsolutions.jump.workbench.ui.plugin.FeatureInstaller;


public class MICPlugIn extends AbstractPlugIn{

    private MICEngine engine = new MICEngine();
    
    public MICPlugIn() {
        // empty constructor
    }

    public String getName() {
        return "Maximum Inscribed Circle"; // set the name / label of plugin for tooltips ...
    } 
    
    private Layer layer;
    String attribute;
    private static String Katmanim = "Select a polygon layer:";
    boolean use_attribute;
  
    public void initialize(PlugInContext context) throws Exception {
        FeatureInstaller featureInstaller = new FeatureInstaller(context.getWorkbenchContext());
   	        featureInstaller.addMainMenuPlugin(this, 
   	        	 new String[] {"Plugins"}, //menu path
   	        	 "Maximum Inscribed Circle (MICGIS)", //name
   	        	 false, //checkbox
	        	 new ImageIcon(this.getClass().getResource("/images/mic.png")), //icon 
   	        	 new MultiEnableCheck().add(context.getCheckFactory().createTaskWindowMustBeActiveCheck()).add(context.getCheckFactory().createAtLeastNLayersMustExistCheck(1)));
   	        	 context.getWorkbenchFrame().getToolBar().addPlugIn(new ImageIcon(getClass().getResource("/images/mict.png")), //icon 
                 this, new MultiEnableCheck().add(context.getCheckFactory().createTaskWindowMustBeActiveCheck()).add(context.getCheckFactory().createAtLeastNLayersMustExistCheck(1)),
                 context.getWorkbenchContext());
    }
    
    public static MultiEnableCheck createEnableCheck(WorkbenchContext workbenchContext) {
        EnableCheckFactory checkFactory = new EnableCheckFactory(workbenchContext);

        return new MultiEnableCheck()
        .add(checkFactory.createWindowWithLayerNamePanelMustBeActiveCheck());
    }

    
    public boolean execute(PlugInContext context) throws Exception {


    	// script for GUI
    	
    	final MultiInputDialog dialog = new MultiInputDialog(context.getWorkbenchFrame(), "Maximum Inscribed Circle Finder", true);
  
    	dialog.setSideBarImage(IconLoader.icon("/images/side.png"));
   	 	dialog.setSideBarDescription("This plugin calculates Maximum Inscribed Circle (MIC) for polygon features by using MICGIS algorithm developed " +
		 		"by Beyhan, Güler and Tağa (2020). Further explanations can be found in the following paper:" + '\n' + '\n' + 
		 		"Beyhan, B., Güler, C. & Tağa, H. (2020) An algorithm for maximum inscribed circle based on Voronoi diagrams and geometrical properties." +
		 		" Journal of Geographical Systems, 22, 391–418."  + '\n' + '\n' + "Source code can be downloaded from https://github.com/burakbeyhan/maximum-inscribed-circle");
   	 
        JComboBox jcb_layer = dialog.addLayerComboBox(Katmanim, context.getCandidateLayer(0), "please select a layer", context.getLayerManager());
        jcb_layer.setPreferredSize(new Dimension(126, 20));

        JCheckBox cikti = dialog.addCheckBox("Use a field for a unique identifier for each PF", false, "for enumaration of PFs select a field");
        List fieldim = getFieldsFromLayerWithoutGeometry(context.getCandidateLayer(0));       
        Object val = fieldim.size() > 0 ? fieldim.iterator().next() : null;
        
        final JComboBox jcb_attribute = dialog.addComboBox("Select a field:", val, fieldim, "please select a field");
        jcb_attribute.setEnabled(false);
        jcb_attribute.setPreferredSize(new Dimension(126, 20));
        
		dialog.addLabel(" ");
		dialog.addLabel("<HTML><EM>"+"Maximum Inscribed Circle Options"+"<HTML><EM>");
   	    dialog.addDoubleField("Core search radius:", 0.75, 11, "search radius");
        JCheckBox MIClayer = dialog.addCheckBox("Add Maximum Inscribed Circles layer", true, "add MICs as a layer");
        JCheckBox MICenter = dialog.addCheckBox("Add Centers of MICs layer", true, "add centers of MICs as a layer");

		dialog.addLabel(" ");
		dialog.addLabel("<HTML><EM>"+"Unit and Precision Specification for Display"+"<HTML><EM>");
		List tips = new ArrayList();
		tips.add("millimeter (mm)");
		tips.add("centimeter (cm)");
		tips.add("decimeter (dm)");
        tips.add("meter (m)");
        tips.add("hectometer (hm)");
        tips.add("kilometer (km)");
  		final JComboBox CB = dialog.addComboBox("Unit of radius and perimeter:", "meter (m)", tips, "select a unit");
  		CB.setPreferredSize(new Dimension(126, 20));

        List cips = new ArrayList();
        cips.add("square mm (mm²)");
		cips.add("square cm (cm²)");
		cips.add("square dm (dm²)");
        cips.add("square m (m²)");
        cips.add("decare (daa)");
        cips.add("hectare (ha)");
        cips.add("square km (km²)");
  		final JComboBox AB = dialog.addComboBox("Unit of surface area:", "square m (m²)", cips, "select a unit");
  		AB.setPreferredSize(new Dimension(126, 20));

  		dialog.addIntegerField("Precision for display:", 15, 11, "specify the precision for displaying results in attribute table");   

  		jcb_layer.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	layer = dialog.getLayer(Katmanim);
                List list = getFieldsFromLayerWithoutGeometry(layer);
                if (list.size() == 0) {
                	jcb_attribute.setModel(new DefaultComboBoxModel(new String[0]));
                	jcb_attribute.setEnabled(true);
                }
                jcb_attribute.setModel(new DefaultComboBoxModel(list.toArray(new String[0])));
            }            
        });

		ActionListener disaver = new ActionListener() {
		      public void actionPerformed(ActionEvent e) {
		        boolean ciktim = dialog.getBoolean("Use a field for a unique identifier for each PF") ;
		        jcb_attribute.setEnabled(ciktim);
		      }
		};
		cikti.addActionListener(disaver);
		        
		ActionListener goruver = new ActionListener() {
		      public void actionPerformed(ActionEvent e) {
		    	boolean micm = dialog.getBoolean("Add Centers of MICs layer") ;
		    	boolean micd = dialog.getBoolean("Add Maximum Inscribed Circles layer") ;
                if (micm && micd) dialog.setSideBarImage(IconLoader.icon("/images/side.png"));
                else if (!micm && micd) dialog.setSideBarImage(IconLoader.icon("/images/sidi.png"));
                else if (micm && !micd) dialog.setSideBarImage(IconLoader.icon("/images/sida.png"));
                else dialog.setSideBarImage(IconLoader.icon("/images/sido.png"));
		      }
		};
		MICenter.addActionListener(goruver);
		MIClayer.addActionListener(goruver);
		
        GUIUtil.centreOnWindow(dialog);
        dialog.setVisible(true);

        if (!dialog.wasOKPressed()) {
            return false;
        }

        getDialogValues(dialog);
        engine.execute(context);

        return true;
    }

    private List getFieldsFromLayerWithoutGeometry(Layer lyr)
    {
      List fields = new ArrayList();
      FeatureSchema schema = lyr.getFeatureCollectionWrapper().getFeatureSchema();
      for (int i = 0; i < schema.getAttributeCount(); i++) {
        if (schema.getAttributeType(i) != AttributeType.GEOMETRY) {
          fields.add(schema.getAttributeName(i));
        }
      }
      return fields;
    }
   
    private void getDialogValues(MultiInputDialog dialog) {
        engine.aktarim1(dialog.getLayer(Katmanim));
        engine.aktarim2(dialog.getText("Select a field:"));
        engine.aktarim3(dialog.getText("Unit of radius and perimeter:"));
        engine.aktarim4(dialog.getBoolean("Add Maximum Inscribed Circles layer"));
        engine.aktarim5(dialog.getBoolean("Add Centers of MICs layer"));
        engine.aktarim6(dialog.getText("Unit of surface area:"));
        engine.aktarim7(dialog.getBoolean("Use a field for a unique identifier for each PF"));
        engine.aktarim8(dialog.getDouble("Core search radius:"));
        engine.aktarim9(dialog.getInteger("Precision for display:"));
    }
    
}

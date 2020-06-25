
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
import com.vividsolutions.jump.workbench.plugin.Extension;
import com.vividsolutions.jump.workbench.plugin.PlugInContext;

/**
 * 
 *  - this class loads the PlugIn into OpenJUMP
 *
 *  - class has to be called "Extension" on the end of classname
 *    to use the PlugIn in OpenJUMP 
 */

public class MICExtension extends Extension{

	/**
	 * calls PlugIn using class method xplugin.initialize() 
	 */
	
    public String getName() {
        return "Maximum Inscribed Circle Finder (Burak Beyhan, Cüneyt Güler and Hidayet Tağa)"; // set the name / label of plugin for tooltips ...
    } 
    
    public String getVersion() {
        return "1.2.0 (2020-06-20)";
    }
    
	public void configure(PlugInContext context) throws Exception{
		new MICPlugIn().initialize(context);
	}
	
}

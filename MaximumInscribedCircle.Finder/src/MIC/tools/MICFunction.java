
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
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import com.vividsolutions.jump.feature.Feature;
import com.vividsolutions.jump.feature.FeatureCollection;
import com.vividsolutions.jump.feature.FeatureDatasetFactory;
import com.vividsolutions.jts.algorithm.Angle;
import com.vividsolutions.jts.algorithm.CGAlgorithms;
import com.vividsolutions.jts.algorithm.MinimumBoundingCircle;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineSegment;
import com.vividsolutions.jts.geom.LinearRing;
import com.vividsolutions.jts.geom.MultiPoint;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.geom.Triangle;
import com.vividsolutions.jts.operation.buffer.BufferOp;
import com.vividsolutions.jts.triangulate.VoronoiDiagramBuilder;


public class MICFunction {

	  // If MIC is delimited by 2 edges and 1 point, or 1 edge and 2 points, following function is used in order to calculate MIC.
	  public static class Apollonius {
	     public static Coordinate getMIC(Coordinate f2, Coordinate s1, Coordinate s2, LineSegment segd, int amod) {
	        Coordinate result = null; 
	        // following part of the function considers the entity configuration of 2 edges and 1 point (001, 010, or 100)
		    if (amod == 1) { // it is assumed that entity configuration is arranged to be presented by the specific case of 100
	    	    LineSegment segc = new LineSegment(s1, s2);
	    	    Coordinate so = segc.midPoint(); // mid point of line c to prevent sequentiality
	    	    Coordinate to = segd.midPoint(); // mid point of line d to prevent sequentiality
	    	    LineSegment segf = new LineSegment(so, to); // line segment between these mid points
	    	    Coordinate fo = segf.midPoint(); // mid point of mid points
	            Coordinate teser = segc.lineIntersection(segd); // computation of the intersection of the lines of infinite extent defined by line entities, segc and segd.
	            double angc = segc.angle();
	            if (teser == null) { // if segc and segd (line entities) are parallel to each other
	            	double rx = segc.distance(segd)/2;
	            	double bxcd = Math.tan(angc); // slope of the mid-line
	            	double axcd = fo.y-bxcd*fo.x;
	    	        double bicd = Math.tan(angc+Math.PI/2); // slope of the line perpendicular to the mid-line
	     	 		double aicd = f2.y - bicd*f2.x;
	        		double xfic = (axcd-aicd)/(bicd-bxcd); // intersection of the mid-line and the line perpendicular to the mid-line and passing through first entity (f2)
	        		double yfic = axcd+bxcd*xfic;
	        		double diko = Math.sqrt((f2.x-xfic)*(f2.x-xfic)+(f2.y-yfic)*(f2.y-yfic)); // distance between the first entity and point of intersection (xfic,yfic)
	        		double yato = Math.sqrt(rx*rx-diko*diko); // distance to the MIC along the mid-line
	            	double dx = yato*Math.cos(angc);
	            	double dy = yato*Math.sin(angc);
	            	Coordinate mic = new Coordinate(xfic+dx, yfic+dy);
	            	Coordinate mci = new Coordinate(xfic-dx, yfic-dy);
	            	if (mic.distance(fo) > mci.distance(fo)) result = mci;
	            	else result = mic;
	            } 
	            else { // if segc and segd are not parallel to each other
	            	double hacim = Angle.angleBetweenOriented(so, teser, to)/2; // calculation of the angle bisector
	            	double cdan = angc+hacim; // angle of the line represented by the angle bisector of the angle formed by the extensions of lines c and d 
	            	if (hacim < 0) cdan = segd.angle()-hacim;;
	            	double dis = (so.distance(teser)+to.distance(teser))/2;
	            	double xis = dis*Math.cos(cdan);
	            	double yis = dis*Math.sin(cdan);
	            	Coordinate tic = new Coordinate(teser.x+xis, teser.y+yis); // a point along the line represented by angle bisector
	            	Coordinate tci = new Coordinate(teser.x-xis, teser.y-yis);
	            	if (tic.distance(fo) > tci.distance(fo)) tic = tci;
	            	double re = segd.distancePerpendicular(tic); // radius of the circle whose center is tic
			   	
	            	// intersection of the line passing through the point entity (f2) and teser with the circle centered around tic
	            	Coordinate yar = null;
	            	double ddx = teser.x-f2.x;
	            	double ddy = teser.y-f2.y;
	            	double aA = ddx*ddx + ddy*ddy;
	            	double bB = 2*(ddx*(f2.x-tic.x) + ddy*(f2.y-tic.y));
	            	double cC = (f2.x-tic.x)*(f2.x-tic.x) + (f2.y-tic.y)*(f2.y-tic.y) - re*re;
	            	double det = bB*bB - 4*aA*cC;
	            	if (det >= 0) {
	            		double t = (-bB+Math.sqrt(det))/(2*aA);
	            		yar = new Coordinate(f2.x+t*ddx, f2.y+t*ddy); // 1st intersection
				   	 	t = (-bB-Math.sqrt(det))/(2*aA);
				   	 	Coordinate sar = new Coordinate(f2.x+t*ddx, f2.y+t*ddy); // 2nd intersection
				   	 	if (yar.distance(teser) < sar.distance(teser)) yar = sar; // working intersection
	            	}

	            	// center of MIC: intersection of the angle bisector with the line having the same slope with the line passing through tic and yar,
	            	if (yar != null ) { // but passing through the point entity (f2)
	            		double acic = Math.atan2(yar.y-tic.y, yar.x-tic.x);
	            		double bmic = Math.tan(acic); // slope of the line that will pass through the first entity (f2)
	            		double amic = f2.y-bmic*f2.x; // intercept of this line
	                	double bcd = Math.tan(cdan); // slope of the angle bisector line (cdan)
	                	double acd = teser.y-bcd*teser.x; // intercept of this line
	            		double xmic = (acd-amic)/(bmic-bcd); // intersection of these lines
	            		double ymic = amic+bmic*xmic;
	            		result = new Coordinate(xmic,ymic);
	            	}
	            }
		    }
		      
	        // following part of the function considers the entity configuration of 1 edge and 2 points (011, 101, or 110)
		    if (amod == 2) { // it is assumed that entity configuration is arranged to be presented by the specific case of 110 
	    	    LineSegment sege = new LineSegment(f2, s2); // creates a line between the first and second point entities
	            Coordinate teser = sege.lineIntersection(segd); // computation of the intersection of the lines of infinite extent defined by segd and sege.
	            Coordinate sto = sege.midPoint(); // if segd and sege are parallel to each other,
	            Coordinate tto = segd.closestPoint(sto); // the closest point to the mid point of sege on segd is used as tto
	            if (teser != null) { // if segd and sege are not parallel to each other,
	            	double fr = Math.sqrt(f2.distance(teser)*s2.distance(teser)); // power of a point theorem is employed to calculate tto
	            	double dx = fr*Math.cos(segd.angle());
	            	double dy = fr*Math.sin(segd.angle());
	            	double xc = teser.x+dx;
	            	double yc = teser.y+dy;
	            	if (!((segd.maxX() > xc && xc > segd.minX()) || (segd.minY() < yc && yc < segd.maxY()))) { // 
	            		xc = teser.x-dx; yc = teser.y-dy;
	            	}
	            	tto = new Coordinate(xc, yc);
	            }
	            result = Triangle.circumcentre(f2, s2, tto);
		    }
		      
		    return result;
	     }
	  }
	  
	  
	  public static double[] getMIC(Polygon poly) { // MIC function with default input parameters ...
	      Envelope penvelope = poly.getEnvelopeInternal(); // for offsetting edges to increase the precision
	      Geometry pcnvxhull = poly.convexHull(); // calculation of convex hull
		  double[] result = MICFunction.getMIC(poly, penvelope, pcnvxhull, 0.75, 2500);
		  return result;
	  }
	  
	  public static double[] getMIC(Polygon poly, double csr, int bpli) { // if csr is determined by the user ...
	      Envelope penvelope = poly.getEnvelopeInternal(); // for offsetting edges to increase the precision
	      Geometry pcnvxhull = poly.convexHull(); // calculation of convex hull
		  double[] result = MICFunction.getMIC(poly, penvelope, pcnvxhull, csr, bpli);
		  return result;
	  }
	  
	  public static double[] getMIC(Polygon poly, Envelope penvelope, Geometry pcnvxhull, double csr, int bpli) { // if envelope and convex hull of the polygon is already calculated ...
	      Polygon geom = (Polygon) poly.clone(); // for processing polygon with holes ...
	      final GeometryFactory gf = new GeometryFactory(); // creating a geometry factory for subsequent calculations
	      int div2 = 2; // used in combination with bpli parameter defining a treshold for the number of edges in splitting the edges of PF concerned into two before considering lpei 
	      double precx = Math.round(penvelope.getMinX())-1; // for offsetting edges to increase the precision
	      double precy = Math.round(penvelope.getMinY())-1; // for offsetting edges to increase the precision

	      // Following part of the script performs two functions: (1) removal of vertices creating collinear segments
	      // and (2) offsetting edges to increase the precision for subsequent calculations.
		  LinearRing geomext = (LinearRing) geom.getExteriorRing(); // consideration of exterior ring for polygons with holes
		  Coordinate[] fas = geomext.getCoordinates();
		  int KEs = fas.length-1;
		  ArrayList<Coordinate> faso = new ArrayList<Coordinate>();
	      faso.add(new Coordinate((fas[0].x-precx), (fas[0].y-precy))); // first vertex of exterior ring is offset and added to the new exterior ring
		  for (int i = 1; i < KEs; i++) {
		     double bddevi = Angle.angleBetween(fas[i-1], fas[i], fas[i+1]); // for removal of vertices creating collinear segments
			 if (bddevi != 0 && bddevi != Math.PI) faso.add(new Coordinate((fas[i].x-precx), (fas[i].y-precy)));
		  }
		  faso.add(new Coordinate((fas[KEs].x-precx), (fas[KEs].y-precy)));
		  Coordinate[] fasi = faso.toArray(new Coordinate[] {});
		  int KSs = fasi.length-1;
		  geomext = gf.createLinearRing(fasi);
		  int doner = geom.getNumInteriorRing(); // number of holes in the polygon
		  LinearRing[] geomina = new LinearRing[doner];
		  for (int ih = 0; ih < doner; ih++) { // consideration of each hole
		   	 fas = geom.getInteriorRingN(ih).getCoordinates();
			 int KEh = fas.length-1;
			 KEs = KEs + KEh;
			 faso = new ArrayList<Coordinate>();
		     faso.add(new Coordinate((fas[0].x-precx), (fas[0].y-precy))); // first vertex of hole considered is offset and added to the new hole concerned
			 for (int hi = 1; hi < KEh; hi++) {
			   	double bdhevi = Angle.angleBetween(fas[hi-1], fas[hi], fas[hi+1]); // for removal of vertices creating collinear segments
			    if (bdhevi != 0 && bdhevi != Math.PI) faso.add(new Coordinate((fas[hi].x-precx), (fas[hi].y-precy)));
			 }
			 faso.add(new Coordinate((fas[KEh].x-precx), (fas[KEh].y-precy)));
			 fasi = faso.toArray(new Coordinate[] {});
			 KSs = KSs + fasi.length-1;
			 LinearRing geomint = gf.createLinearRing(fasi);
			 geomina[ih] = geomint;
		  }
		  geom = gf.createPolygon(geomext, geomina);

		  // For consideration of polygons with holes in the subsequent parts of the main script following array list is created,
	      doner = geom.getNumInteriorRing(); 
	      ArrayList<Integer> geomit = new ArrayList<Integer>();
	      geomit.add(0);
	      if(doner > 0) {
	         geomit.add(geom.getExteriorRing().getCoordinates().length);
		     for (int ih = 0; ih < doner; ih++) geomit.add(geomit.get(ih)+geom.getInteriorRingN(ih).getCoordinates().length);
	      }

	      double P = geom.getLength(); // calculation of perimeter of polygon
	      double A = geom.getArea(); // calculation of area of polygon

	      
	      // Calculation and creation of maximum inscribed circles and writing MIC parameters to attribute tables.
		  List<Point> points = new ArrayList<Point>();
		  Coordinate[] as = geom.getCoordinates();
		  int aslen = as.length; // total number of nodes in the polygon
	      geomit.add(aslen);
		  double edgbl = 2.0000001; // for creation of mid points
	      int NEs = aslen-1; // total number of edges in polygon
	  	  double AL = P/NEs; // calculation of average edge length
	      Coordinate bbkng; 
	      double bbimp;
	      if (NEs > bpli) div2 = 1; // if the number of edges is more than bpli, the edges of the PF is not automatically divided into two parts.
		    						// this fastens the calculation of MICs for big polygons having huge number of edges
	      
		  LinearRing[] geomlin = new LinearRing[doner+1];
		  geomlin[0] = (LinearRing) geom.getExteriorRing();
		  for (int ih = 0; ih < doner; ih++) {
			 LinearRing geomint = (LinearRing) geom.getInteriorRingN(ih);
			 geomlin[ih+1] = geomint;
		  }

		  if (NEs > 2 ) { // for extraction of vertices of polygon and population of input points that will be used in the creation of Voronoi diagram
			 for (int lr = 0; lr < doner+1; lr++) { // for consideration of polygon with holes in this process ...
			    LinearRing geomln = geomlin[lr];
			    Coordinate[] kas = geomln.getCoordinates();
		        int NEr = kas.length-1;
		        for (int i = 0; i < NEr; i++) { // iteration over nodes and edges  
		    		com.vividsolutions.jts.geom.Point pt = gf.createPoint(kas[i]);
		    		points.add(pt);
		    		Coordinate p1 = kas[i];
		    		Coordinate p2 = kas[i+1];
		    		double ELx = (p2.x-p1.x);
		    		double ELy = (p2.y-p1.y);
		    		double EL = Math.sqrt(ELx*ELx+ELy*ELy);
		    		if (div2 == 2) points.add(gf.createPoint(new Coordinate(ELx/edgbl+p1.x, ELy/edgbl+p1.y)));
		    		double bolen = Math.round(EL/AL); // calculation of bolen as a parameter for densification of input points
		    		double lpei = (EL/bolen)/div2; // calculation of length of point enrichment interval (lpei)
		    		double xlpei = lpei*(ELx/EL);
		    		double ylpei = lpei*(ELy/EL);
		  	    	for (int j = 1; j < bolen; j++) { // for enrichment of points by using lpei 
		    			points.add(gf.createPoint(new Coordinate(p1.x+xlpei*j, p1.y+ylpei*j)));
		    			if (div2 == 2) points.add(gf.createPoint(new Coordinate(ELx/edgbl+p1.x+xlpei*j, ELy/edgbl+p1.y+ylpei*j)));
		  	       } 
		    	}
			 }
		  }	
		    	
		  MultiPoint allpoints = gf.createMultiPoint(points.toArray(new Point[points.size()])); // creation of multipoints for the construction of Voronoi diagram
		  VoronoiDiagramBuilder voronoiBuilder = new VoronoiDiagramBuilder(); // construction of Voronoi diagram
		  voronoiBuilder.setSites(allpoints);
		  voronoiBuilder.setTolerance(0);
		  Geometry gv = voronoiBuilder.getDiagram(gf); // Voronoi diagram constructed
	      ArrayList<Geometry> gvlist = new ArrayList<Geometry>();
	      for (int i = 0; i < gv.getNumGeometries(); i++) gvlist.add(gv.getGeometryN(i));
	      FeatureCollection gvcol = FeatureDatasetFactory.createFromGeometry(gvlist);

	      // Extraction of the approximate Medial Axis (MA) of the polygon via selection of the vertices of Voronoi diagram inside the polygon and
		  ArrayList<Geometry> points1 = new ArrayList<Geometry>(); // enrichment of them by adding mid points
		  for (Iterator ic = gvcol.iterator(); ic.hasNext();) {
		     Feature ice = (Feature) ic.next();
		     Geometry ici = ice.getGeometry();
		     Coordinate[] vas = ici.getCoordinates();
		     for (int i = 0; i < vas.length-1; i++) {
		        com.vividsolutions.jts.geom.Point pta = gf.createPoint(vas[i]);
		        com.vividsolutions.jts.geom.Point ptb = gf.createPoint(new Coordinate((vas[i].x+vas[i+1].x)/2,(vas[i].y+vas[i+1].y)/2)); // creating mid point between the subsequent vertices
		        if (pta.within(geom)) points1.add(pta); // adding vertex as a point to points1 if it is inside the polygon
		        if (ptb.within(geom)) points1.add(ptb); // adding mid point as a point to points1 if it is inside the polygon
		     }
		  }
		        
	      Point geoc = geom.getCentroid(); // centroid of the polygon
	      if(geoc.within(geom)) points1.add(geoc); // adding the centroid as a point to points1 if it is inside the polygon 
		        
		  HashSet hs = new HashSet();
		  hs.addAll(points1);
		  points1.clear();
		  points1.addAll(hs);
		  int cndts = points1.size();
		        
		  double extpole = Math.sqrt(csr); // extra parameterization of csr used in core search
		  
		  // Treatment of very irregular polygons for which Voronoi diagram could not result in adequate number of points in points1
		  double marator = pcnvxhull.getArea(); // using the area of convex hull as an additional parameter for the identification of these irregular polygons 
		  if(cndts <= 3 || marator/A > 8) { // selection of these polygons on the base of the number of points in points1 and the area of convex hull of the polygon
		     extpole = 0.5; // for these irregular polygons extpole is also reduced to cover more inner points as alternative cores ...
		     Point geoi = geom.getInteriorPoint(); // interior point of the polygon
		     if(!geoi.equals(geoc)) points1.add(geoi); // adding interior point to points1 if it is different from the centroid of the polygon
	 		 for (int lr = 0; lr < doner+1; lr++) { // consideration for polygon with holes ... 
				LinearRing geomln = geomlin[lr];
				Coordinate[] kas = geomln.getCoordinates();
				int NEr = kas.length-1;
			    for (int i = 0; i < NEr; i++) {
				   Coordinate p1 = kas[i];
				   Coordinate p2 = kas[i+1];
				   int has = i+2;
			       if (has >= NEr) has = has - NEr;
				   Coordinate p3 = kas[has];
				   double EL = p1.distance(p2);
				   double EK = p2.distance(p3);
				   double divo = edgbl;
				   double pivo = (EK/EL)*divo;
				   if(EL > EK) {pivo = edgbl; divo = (EL/EK)*pivo;}
				   Coordinate arpnt = new Coordinate(p2.x-(p2.x-p1.x)/divo,p2.y-(p2.y-p1.y)/divo);
				   Coordinate arpns = new Coordinate((p3.x-p2.x)/pivo+p2.x,(p3.y-p2.y)/pivo+p2.y);
				   Coordinate arp = new Coordinate((arpns.x-arpnt.x)/edgbl+arpnt.x,(arpns.y-arpnt.y)/edgbl+arpnt.y);
				   com.vividsolutions.jts.geom.Point ptra = gf.createPoint(arp);
				   if(ptra.within(geom)) points1.add(ptra); // if it is inside the polygon 
				}
			 }
		     HashSet hst = new HashSet();
			 hst.addAll(points1);
			 points1.clear();
			 points1.addAll(hst);
		  }
		  cndts = points1.size(); // final number of points in points1 ...
		        
		  // Calculation of the first core and the first candidate for the center of MIC together with the closest first three edges to this center
		  MinimumBoundingCircle mbcu = new MinimumBoundingCircle(geom);
		  double distemp = mbcu.getRadius()*2;
		  double[] mindist1 = new double[cndts]; // for storage of the minimum distance between the Voronoi vertex and the edges of the polygon
		  int[] segment1 = new int[cndts]; // for storage of id of the first node of the edge closest to the voronoi vertex concerned
		  double[] mindist2 = new double[cndts];
		  int[] segment2 = new int[cndts];
		  double[] mindist3 = new double[cndts];
		  int[] segment3 = new int[cndts];
		  for (int j = 0; j < cndts; j++) {
			 double mes1 = distemp;
			 double mes2 = distemp;
			 double mes3 = distemp;
			 int ed1 = -1;
			 int ed2 = -1;
			 int ed3 = -1;
		     Coordinate p0 = points1.get(j).getCoordinate();
			 for (int lr = 0; lr < doner+1; lr++) { // consideration for polygon with holes ...
				int eki = geomit.get(lr);
				LinearRing geomln = geomlin[lr];
				Coordinate[] kas = geomln.getCoordinates();
				int NEr = kas.length-1;
				int aod = NEr-1;
			    for (int i = 0; i < NEr; i++) {
			       int adi = i;
		           double dist = CGAlgorithms.distancePointLine(p0, kas[i], kas[i+1]);
		           double disl = CGAlgorithms.distancePointLine(p0, kas[aod], kas[i]);
		           if (disl == dist) adi = aod;
		           if (dist <= mes1) {
		           	 mes3 = mes2; mes2 = mes1; mes1 = dist; 
		           	 ed3 = ed2; ed2 = ed1; ed1 = adi + eki;
		           }
		           else if (dist <= mes2 && i+eki != ed1) {
		           	 mes3 = mes2; mes2 = dist; 
		           	 ed3 = ed2; ed2 = adi + eki;
		           }
		           else if (dist <= mes3 && i+eki != ed1 &&  i+eki != ed2) {
		           	 mes3 = dist; 
		           	 ed3 = adi + eki;
		           }
				   aod = i;
		        }
			 }
		     mindist1[j] = mes1;
		     segment1[j] = ed1;
		     mindist2[j] = mes2;
		     segment2[j] = ed2;
		     mindist3[j] = mes3;
		     segment3[j] = ed3;
		  }

		  // Calculation of the radius of the first MIC candidate and id of the first node of the edge tangent to the respective MIC
		  double mr = 0;
		  int index = 0;
		  int edi1 = 0;
		  for (int j = 0; j < cndts; j++) {
		     if (mindist1[j] > mr) {
		        mr = mindist1[j]; // radius (mr) of the first MIC candidate
		        index = j; // point id (n1) of the MIC concerned in points1
		        edi1 = segment1[j]; // id of the first vertex of the edge (a) tangent to the MIC
		     }
		  }
		        
		  Geometry ilkcor = points1.get(index); // the first MIC candidate (n1)
		  double micCircleRad = mr;
		  double micCirclecheck = 0;
		  double micX = ilkcor.getCoordinate().x;
		  double micY = ilkcor.getCoordinate().y;
		  double micXcheck = micX;
		  double micYcheck = micY;
	  	  int ncore = 1;
	  		
	  	  // Finding the closest point (n2) in points1 to n1 among the points having a distance to edge (a) more than the radius (mr) of first MIC candidate
		  Coordinate p0 = ilkcor.getCoordinate();
		  double dists = distemp;
		  int indes = 0;
		  for (int j = 0; j < cndts; j++) {
		     Coordinate p1 = points1.get(j).getCoordinate();
		     double dist = p0.distance(p1); // distance between n1 and n2
		     double disl = CGAlgorithms.distancePointLine(p1, as[edi1], as[edi1+1]); // distance between n2 and edge (a)
		     if(dist < dists && j != index && mr < disl) {
		        dists = dist;
		        indes = j;
		     }
		  }
		  double disto = mindist1[indes]; // the distance between n2 and the edge closest to it
		        
		  // Searching for the potential cores that may cover the center of MIC
		  Comparator<double[]> doucomp = new Comparator<double[]>() { // creating a comparator for sorting 
			 public int compare(double[] a, double[] b) { // a double array according to the first column 
			    return Double.compare(b[0], a[0]); // in descending order
			 }
		  };
		  ArrayList canodam = new ArrayList(); // creating a list (canodam) to store information about the points for canodak
		  for (int j = 0; j < cndts; j++) { // all the points in points1 are checked against whether the distance between them and the edge closest to them 
			 if (mindist1[j] > disto*extpole) canodam.add(j); // is more than the multiplication of the distance (disto) between n2 and the edge closest to it 
		  } // with the square root of csr (here extpole)
		  double[][] canodak = new double[canodam.size()][4]; // putting candom into canodak together with the other required information
		  for (int j = 0; j < canodak.length; j++) {
			 int jj = (Integer) canodam.get(j);
			 canodak[j][0] = mindist1[jj];
		   	 canodak[j][1] = (double) segment1[jj];
		   	 canodak[j][2] = (double) jj;
		  }
		  Arrays.sort(canodak, doucomp); // sorting canodak according to the first column in descending order
		  int ceper = 0;
		  for (int j = 0; j < canodak.length; j++) { // canodak is further processed in order to find the centers of the cores;
			 Geometry pc1 = points1.get((int) canodak[j][2]); // creating the initial point from the information for this point in canodak
			 for (int k = 0; k < canodak.length; k++) { 
			 	if (j != k && canodak[j][3] == 0 && canodak[k][3] == 0) {
			 		Geometry pc2 = points1.get((int) canodak[k][2]); // creating the subsequent point from the information for the point concerned in canodak
			 		if (pc2.isWithinDistance(pc1, mr*csr)) { // those points that are included within the circle drawn around n1 (first pc1) with a radius
			 			canodak[k][3] = 1; // defined according to the multiplication of mr with csr, and subsequently other core candidates (subsequent pc1s)
			 			ceper++; // are identified and marked with 1 in the forth column of canodak (1 indicates that the point concerned is not a possible core)
			 		}
			 	}
			 }
		  }
		  double[][] cod = new double[canodak.length-ceper][3]; // those points not marked with 1 are stored in a new array (cod) representing centers of possible cores
		  ceper = 0;
		  for (int j = 0; j < canodak.length; j++) {
			 if (canodak[j][3] == 0) { // those points marked with 1 are not considered as a possible core
				cod[ceper][0] = canodak[j][0];
				cod[ceper][1] = canodak[j][1];
				cod[ceper][2] = canodak[j][2];
		    	ceper++; 
			 }
		  }

		  // Consideration of all core alternatives in the delimitation of MIC
		  ncore = ceper; // total number of cores 
		  for (int colp = 0; colp < ncore; colp++) { 
			 double maxdist = cod[colp][0]; // getting information about the core
			 edi1 = (int) cod[colp][1]; // starting point id for the first edge delimiting approximate MIC
			 index = (int) cod[colp][2]; // point id for the approximate center of MIC located along the medial axis
			 ilkcor = points1.get(index); // the approximate center of MIC
			 Coordinate mic1 = ilkcor.getCoordinate();

			 ArrayList edi1alt = new ArrayList(); // creation of a list for limited edi1 alternatives
			 int edi2to1;
			 int ed2i = segment2[index]; // starting point id for the second closest edge to the approximate center of MIC
			 int ed3i = segment3[index]; // starting point id for the third closest edge to the approximate center of MIC
			 
			 // Consideration of limited alternatives for the first edge (edi1)
			 for (int felp = 0; felp < 3; felp++) {
		        Coordinate n1 = ilkcor.getCoordinate();
		        
	            int NEh = NEs;
	            for (int aos = 1; aos < geomit.size(); aos++) {
	            	if(edi1 >= geomit.get(aos-1) && edi1 <= geomit.get(aos)) NEh = geomit.get(aos)-1;
	            }
		        
	            int edis = edi1+1; // ending point id for the first edge delimiting approximate MIC
		    	if (edis >= NEh) edis = edis - NEh;
		        Coordinate n3 = as[edi1]; // starting point of the first edge
		        Coordinate n4 = as[edis]; // ending point of the first edge
	    	    LineSegment sega = new LineSegment(n3, n4); // construction of the first edge (a) delimiting MIC
		        double b22 = sega.angle(); // computation of the angle that the vector defined by this edge makes with the X-axis 

		        int fedge = 0; // if the distance between n1 and the ending node of the first edge is equal to the distance between n1 and the edge concerned,
		        if (n4.equals(sega.closestPoint(n1))) fedge = 1;  // it is marked with a possibility of being touched by the MIC at the ending node concerned.
		        // in other words, if the closest point to n1 on the first edge is the ending node of the edge concerned, it is marked with 1 for Apollonius function.
	            
	            // Consideration of alternatives for edi2 (second edge)
			    ArrayList<double[]> edi2alt = new ArrayList(); // creation of an array list to store the alternatives for edi2
		        double bi2 = Math.tan(b22 + Math.PI/2); // slope of the line perpendicular to the first edge (line a)
	 	 		double ako = n1.y - bi2*n1.x; // if a point located along this line and passing through n1 is positioned at the reverse side of line a at a certain
		        double cbck = Math.cos(b22 + Math.PI/2); // distance from n1, the closest edge to the point concerned can be designated as the 2nd edge (Type 2).
		        double deltx = cbck*(mr/Math.PI); // the distance between this point and n1 is defined by dividing mr by Ï€
		        Coordinate mp1 = new Coordinate(n1.x-deltx, ako+bi2*(n1.x-deltx)); // the first alternative for the point concerned - default alternative
		        Coordinate mp2 = new Coordinate(n1.x+deltx, ako+bi2*(n1.x+deltx)); // the second alternative for the point concerned
		        if (sega.distance(mp2) > sega.distance(mp1)) mp1 = mp2; // if mp2 is further from the first edge, it is designated as the valid alternative.

		        int edis2 = -1; // finding the closest edge (edis2) to mp2 as an alternative for the second edge
		        double dist2 = distemp;
				for (int lr = 0; lr < doner+1; lr++) { // consideration for polygon with holes ...
					int eki = geomit.get(lr);
					LinearRing geomln = geomlin[lr];
					Coordinate[] kas = geomln.getCoordinates();
					int NEr = kas.length-1;
					for (int im = 0; im < NEr; im++) {
			            double distm = CGAlgorithms.distancePointLine(mp1, kas[im], kas[im+1]);
			            if(distm < dist2 && im+eki != edi1) {
			            	dist2 = distm; 
			            	edis2 = im+eki;
			            }
			      	}
				}
			    if (edis2 != -1 && edis2 != edi1) { // if there is an edis2 and it is different from edi1, adding it into edi2alt;
			    	double[] suor = new double[3];
		            for (int aos = 0; aos < geomit.size()-1; aos++) { // construction of edges are based on the starting point id and ending point id. 
		            	int NEp = geomit.get(aos); // those edges whose starting point id is 0 or the starting node of an inner hole does not obey the general 
		            	int NEd = geomit.get(aos+1)-1; // rule of clockwise movement in the construction of the edges. thus all the alternatives for 
		            	if (edis2 == NEp) { // the second edge is checked against to this and necessary correction is done in their starting point id.
		            		edis2 = NEd-1; // this is also considered at the very beginning in segment1, 2 and 3 that are heavily used in the script.
		            	}
		            }
		          	suor[0] = edis2; // alternative for edi2
		          	suor[1] = mp1.x; // x coordinate for n2
		          	suor[2] = mp1.y; // y coordinate for n2
		          	edi2alt.add(suor);
			    }
			      
		        if (ed2i != -1 && ed2i != edis2 && ed2i != edi1) { // if ed2i is different from edis2 and edi1, it is also added into edi2alt;
		        	double[] suor = new double[3];
			        suor[0] = ed2i;
			        suor[1] = (((as[ed2i].x+as[ed2i+1].x)/2)+n1.x)/2;
			        suor[2] = (((as[ed2i].y+as[ed2i+1].y)/2)+n1.y)/2;
			        edi2alt.add(suor);
		        }
		          
		        // consideration for limited alternatives for edi1 on the base of ed2i and edis2
		        if (CGAlgorithms.distancePointLine(n1, as[ed2i], as[ed2i+1]) < CGAlgorithms.distancePointLine(n1, as[edis2], as[edis2+1])) edi1alt.add(ed2i);
		        else edi1alt.add(edis2); // if edis2 is closer to n1 compared with ed2i, it is added into the alternatives for the first edge (edi1alt)
		        edi1alt.add(ed3i); // adding ed3i as an alternative for the first edge into edi1alt
		          
		        int edim3 = -1; int edim4 = -1;
		        int edis3 = -1; int edis4 = -1;
		        double dist3 = maxdist;
		        double dist4 = distemp; // the closest point to n1 in points1 can be declared as n2; alternatively,
		        for (int j = 0; j < cndts; j++) { // among the points located within the first MIC candidate and along the MA, the point furthest from
		            Coordinate pj = points1.get(j).getCoordinate(); // the first edge with a distance of at least mr can also be considered as n2.
		            double dist = n1.distance(pj); // distance between this point and n1
		            double disl = sega.distancePerpendicular(pj); // perpendicular distance between this point and the first edge
		            if(dist < dist4 && j != index) { // finding the closest point to n1 among the points located along the MA (points1)
		            	dist4 = dist; edim4 = j; // id of the point concerned
		            }
		            if(dist3 < disl && j != index && dist < mr) { // finding the point furthest from the first edge but within the radius of circle drawn around n1
		            	dist3 = disl; edim3 = j; // id of the point concerned
		            }
		        }
		        if(edim3 != -1) edis3 = segment1[edim3]; // starting point id of the edge closest to the point represented by edim3 along the MA 
		        if(edis3 == edi1) edis3 = segment2[edim3]; // if this point is identical to the one of the first edge, second closest edge to edim3 is considered
		        if(edim4 != -1) edis4 = segment1[edim4]; // starting point id of the edge closest to the point represented by edim4 along the MA 
		        if(edis4 == edi1) edis4 = segment2[edim4]; // if this point is identical to the one of the first edge, second closest edge to edim4 is considered
			    if (edim3 != -1 && edis3 != edis2 && edis3 != ed2i) { // if edim3 is different from the other alternatives for the second edge, 
			    	double[] suor = new double[3]; // it is added into edi2alt (list of alternatives for the second edge).
			    	suor[0] = edis3;
			    	suor[1] = points1.get(edim3).getCoordinate().x;
			    	suor[2] = points1.get(edim3).getCoordinate().y;
			    	edi2alt.add(suor);
		        }
				if (edim4 != -1 && edis4 != edis2 && edis4 != ed2i && edis4 != edis3) { // if edim4 is different from the other alternatives for the second edge, 
					double[] suor = new double[3]; // it is added into edi2alt.
					suor[0] = edis4;
					suor[1] = points1.get(edim4).getCoordinate().x;
					suor[2] = points1.get(edim4).getCoordinate().y;
					edi2alt.add(suor);
				}

				// Iterating over the alternatives for the second edge
			 	for (int selp = 0; selp < edi2alt.size(); selp++) {
			 		double[] fuor = edi2alt.get(selp);
			 		if (edi1 == (int) fuor[0]) { // extra control for prevention of the identicality between the first and second edges
			 			selp++;
			 			if (selp < edi2alt.size()) fuor = edi2alt.get(selp);
			 			else break;
			 		}
			 		
		    		Coordinate n2 = new Coordinate(fuor[1], fuor[2]);
	        	    LineSegment segb = new LineSegment(n1, n2); // construction of the line (b) connecting approximate center of MIC and n2 
	                double b11 = segb.angle(); // computation of the angle that the vector defined by this line makes with the X-axis 
	                double b1 = Math.tan(b11);
	                double a1 = n1.y - b1*n1.x;

	                int edi2 = (int) fuor[0]; // starting point id for the second edge delimiting approximate MIC

		            NEh = NEs;
		            for (int aos = 1; aos < geomit.size(); aos++) {
		              	if(edi2 >= geomit.get(aos-1) && edi2 <= geomit.get(aos)) NEh = geomit.get(aos)-1;
		            }
		            
	                edis = edi2+1; // ending point id for the second edge delimiting approximate MIC
			    	if (edis >= NEh) edis = edis - NEh;
	                int edit = edi2+2;
			    	if (edit >= NEh) edit = edit - NEh;
	    	        Coordinate n5 = as[edi2]; // starting point of the second edge
	    	        Coordinate n6 = as[edis]; // ending point of the second edge
	        	    LineSegment segc = new LineSegment(n5, n6); // construction of the second edge (c) delimiting MIC
	                double b33 = segc.angle(); // computation of the angle that the vector defined by this edge makes with the X-axis 

	    	        int sedge = 0; // if the closest point to n1 on the second edge is the ending node of the edge concerned, 
	    	        if (n6.equals(segc.closestPoint(n1))) sedge = 1; // it is marked with 1 for Apollonius function
		        
	                ArrayList<int[]> edglist = new ArrayList(); // construction of a list to store the possible configurations of the first and second edges
	                int[] ador = new int[2]; // in relation to whether approximate MIC is tangent to these edges or touching their ending nodes
	                ador[0] = 0; // first edge 
	                ador[1] = 0; // second edge
	                edglist.add(ador); // 0 stands for the fact that the edge concerned is tangent to the MIC and 1 indicates that MIC is touching the ending node of the edge concerned

	                if (fedge != 0 || sedge != 0) { // if the closest point to n1 on the the first or second edge is the ending node of the edge concerned,
	                	int[] edor = new int[2]; // this configuration is added into edglist for consideration of Apollonius function
	                	edor[0] = fedge;
	                	edor[1] = sedge;
	                	edglist.add(edor);
	                }
	            
	                // if the line passing through n2 and perpendicular to the first edge is not crossing the first edge or the line
	                // passing through n1 and perpendicular to the second edge is not crossing the second edge,
	                Coordinate btest = sega.closestPoint(n2); // the edges concerned are again marked with 1.
	                int iffedge = 0; // in other words, if the closest point to n2 on the first edge is one of the ending nodes of the edge concerned,
	                if (n3.equals(btest) || n4.equals(btest)) iffedge = 1; // it is marked with 1.
	                Coordinate ctest = segc.closestPoint(n1);
	                int ifsedge = 0; // and if the closest point to n1 on the second edge is one of the ending nodes of the edge concerned,
	                if (n5.equals(ctest) || n6.equals(ctest)) ifsedge = 1; // it is marked with 1.
	                if (iffedge != fedge || ifsedge != sedge) { // based on the calculation done above, if there is a different configuration from the first one,
		            	int[] edor = new int[2]; // it is also added into edglist for consideration of Apollonius function
		            	edor[0] = iffedge;
		            	edor[1] = ifsedge;
		 		        edglist.add(edor);
	                }

	                // calculations for the second approximation for the center of MIC based on the fact that if MIC is delimited by the first and second edges, the center
	                // of MIC should lie at the intersection of segb and the angle bisector of the angle formed by the extensions of first and second edges (sega and segc).
	                Coordinate keser = sega.lineIntersection(segc); // computation of the intersection of the lines of infinite extent defined by sega and segc.
	                double hacim = 0;
	                double xi = (n3.x+n5.x)/2; 
	                double yi = (n3.y+n5.y)/2;
	                if (keser != null) hacim = Angle.angleBetweenOriented(sega.midPoint(), keser, segc.midPoint())/2; // the angle bisector of the lines of sega and segc
	                if (!(Math.abs(hacim) == 0)) { // if the first and second edges (sega and segc) are not parallel to each other
	                	xi = keser.x;
	                    yi = keser.y;
	                }
	                
	                double an1 = b22 + hacim; // computation of the angle that the vector defined by the angle bisector of the lines of sega and segc makes with the X-axis 
	                if (hacim < 0) an1 = b33 - hacim; 
	                double bck = Math.tan(an1); // slope of the angle bisector of the lines of sega and segc
	                double ack = yi - bck*xi; // intercept of the line represented by the angle bisector of the lines of sega and segc
	                double x1 = (ack-a1)/(b1-bck); // the intersection of segb and the angle bisector of the angle formed by the extensions of sega and segc. 
	                double y1 = a1 + b1*x1;
	                if (b1 == bck) { // if the angle bisector is parallel to segb connecting n1 and n2,
	                	x1 = segb.midPoint().x; // mid point of segb can be taken as the center of the MIC
	                	y1 = segb.midPoint().y;
	                }  
	                
	                Coordinate mic2 = new Coordinate(x1, y1); // coordinate of the second approximation for the center of MIC
	                com.vividsolutions.jts.geom.Point pic2 = gf.createPoint(mic2); // the second approximation for the center of MIC
	                double rad2 = sega.distance(mic2); // radius of the MIC concerned
	                if (!pic2.within(geom)) rad2 = 0;
	                else { // updating micCircleRad, micX and micY if the resulting MIC is better than the existing MIC
	        			bbkng = mic2; // by re-finding the distance between the closest edge to mic2 and mic2
	        			bbimp = rad2;
	        			for (int lr = 0; lr < doner+1; lr++) { // consideration for polygon with holes ...
	        				LinearRing geomln = geomlin[lr];
	        				Coordinate[] kas = geomln.getCoordinates();
	        				int NEr = kas.length-1;
	        				for (int ip = 0; ip < NEr; ip++) {
	        					double fonlem = CGAlgorithms.distancePointLine(bbkng, kas[ip], kas[ip+1]);
	        					if(fonlem <= bbimp) bbimp = fonlem; 
	        				}
	        			}
	        			rad2 = bbimp;
	        			if (bbimp >= micCircleRad) {
	        				micCircleRad = rad2;
	        				micX = mic2.x;
	        				micY = mic2.y;
	        			}
	        		}
	        		
	                // if the line perpendicular to segc and passing through mic2 does not cross segc, the edge crossed by this line  
	                // and contiguous to segc is also considered as an alternative for the second edge (segc).
	                Coordinate nc6 = as[edit]; // for this purpose, firstly edge adjacent to segc should be selected  
	        	    LineSegment segj = new LineSegment(n6, nc6);
	        	    if (segj.distance(mic2) == segj.distancePerpendicular(mic2)) { // if distance between mic2 and segj is equal to perpendicular distance of segj to mic2,
	                	int jest = 0; // it means that the line passing through mic2 and perpendicular to segj is crossing segj that can be considered as a segc alternative.
	                	for (int oci = 0; oci < edi2alt.size(); oci++) { // however, other alternatives for the second edge should be reviewed in order to control whether  
	                		double[] juor = edi2alt.get(oci); // or not new alternative is already included in edi2alt.
	                		int jod2 = (int) juor[0];
	                		if (jod2 == edis) jest = 1; // starting node of the new alternative edge is edis ...
	                	}
	                	if (jest == 0) { // if the new alternative is not included in edi2alt,
	                		double[] suor = new double[3]; // it is added as an alternative for the second edge into edi2alt.
	                		suor[0] = edis;
	                		suor[1] = fuor[1];
	                		suor[2] = fuor[2];
	 				   		edi2alt.add(suor);
	                	}
	                }
	          
	                // Consideration of alternatives for the third edge (edi3) delimiting MIC  
	                ArrayList<Integer> edi3opt = new ArrayList(); // creation of an array list to store the options for edi3
	        	    Coordinate mo1 = sega.closestPoint(mic2); // the first option for the third edge is to select the edge on the base of the first two edges 
	        	    Coordinate mo2 = segc.closestPoint(mic2); // for this purpose, firstly, points at which sega and segc tangent to MIC are calculated.
	                deltx = Math.cos(an1)*(mr/Math.PI); // subsequently, at a distance (defined by dividing mr by Ï€) from mic2, 
	                mp1 = new Coordinate(x1-deltx, ack+bck*(x1-deltx)); // two points (mp1 and mp2) located along the angle bisector of sega and segc are created
	                mp2 = new Coordinate(x1+deltx, ack+bck*(x1+deltx)); // as alternatives that can be used for selection of the alternative edge(s) for edi3 such
	                double m1dst = CGAlgorithms.distancePointLine(mp1, mo1, mo2); // that the closest edge to the point which is further from the circle beam
	                double m2dst = CGAlgorithms.distancePointLine(mp2, mo1, mo2); // connecting the 1st edge and 2nd edge can be used as an alternative for edi3
		        
	                // finding the third closest edge to mic2 and alternatives to this edge on the base of mp1 and mp2
	                int edo3 = -1;
	                int eda3 = -1;
	                int edu3 = -1;
	                double distl = distemp;
	                double dastl = distemp;
	                double dustl = distemp;
	                Coordinate p1;
	                Coordinate p2;
	                for (int lr = 0; lr < doner+1; lr++) { // consideration for polygon with holes ...
	                	int eki = geomit.get(lr);
	                	LinearRing geomln = geomlin[lr];
	                	Coordinate[] kas = geomln.getCoordinates();
	                	int NEr = kas.length-1;
	                	for (int i = 0; i < NEr; i++) {
	                		p1 = kas[i];
	                		p2 = kas[i+1];
	                		double dist = CGAlgorithms.distancePointLine(mic2, p1, p2); // the distance between mic2 and the edge concerned
	                		double dast = CGAlgorithms.distancePointLine(mp1, p1, p2); // the distance between mp1 and the edge concerned
	                		double dust = CGAlgorithms.distancePointLine(mp2, p1, p2); // the distance between mp2 and the edge concerned
	                		if ((dist < distl && i+eki != edi1 && i+eki != edi2) ) { // the closest edge to mic2 can be designated as edi3 
	                			distl = dist; // provided that it is not identical to edi1 and edi2
	                			edo3 = i+eki;
	                		}
	                		if ((dast < dastl && i+eki != edi1 && i+eki != edi2) ) { // the closest edge to mp1 can be considered as an alternative for edi3 
	                			dastl = dast; // provided that it is not identical to edi1 and edi2
	                			eda3 = i+eki;
	                		}
	                		if ((dust < dustl && i+eki != edi1 && i+eki != edi2) ) { // the closest edge to mp2 can also be considered as an alternative for edi3 
	                			dustl = dust; // provided that it is not identical to edi1 and edi2
	                			edu3 = i+eki;
	                		}
	                	}
	                }

	                int mfod = -1; // based on the distances of mp1 and mp2 to the circle beam connecting mo1 and mo2, eda3 or edu3 can be included 
	                int mfid = -1; // in the list of options for edi3 (edi3opt) as mfod. however, if the angle formed by mo1 (point of tangency of the first edge
	                int mwod = -1; // to MIC), mic2 and mo2 (point of tangency of the second edge) is above a certain level (145), both eda3 or edu3 can be used 
	                if (m2dst > m1dst) mfod = edu3; // as alternatives for edi3.
	                else mfod = eda3;
	                double fodvar = Math.toDegrees(Angle.angleBetween(mo1, mic2, mo2));
	                if (Math.abs(fodvar) > 145) { 
	                	if (m2dst > m1dst) mfid = eda3;
	                	else mfid = edu3;
	                }

	                edi3opt.add(mfod); // adding mfod into the options for edi3
	                if (mfid != mfod) edi3opt.add(mfid); // adding mfid into the options for edi3 if it is different from mfod
	                if (edo3 != mfod && edo3 != mfid) edi3opt.add(edo3); // adding edo3 into the options for edi3 if it is different from mfod and mfid
	                if (ed3i != mfod && ed3i != mfid && ed3i != edo3) edi3opt.add(ed3i); // adding ed3i into the options for edi3 if it is different from other options

	                for (int aos = 0; aos < geomit.size()-1; aos++) { // construction of edges are based on the starting point id and ending point id. 
	                	int NEp = geomit.get(aos); // if MIC is not tangent to edge but passing through starting or ending nodes, those edges whose starting point 
	                	int NEd = geomit.get(aos+1)-1; // id is 0 or the starting node of an inner hole does not obey the general rule of clockwise movement in the
	                	if (mfod == NEp || mfid == NEp || edo3 == NEp || ed3i == NEp) { // construction of the edges. thus all the alternatives for the third edge
	                		mwod = NEd-1; // is checked against to this and necessary correction is done in their starting point id.
	                		if (mwod != mfid && mwod != edo3 && mwod != ed3i && mwod != mfod) edi3opt.add(mwod); // this is also important for ed1i and ed2i.
	                	} // thus, this is considered at the very beginning in segment1, segment2 and segment3 that are heavily used in the script.
	                }
				  
	                ArrayList<Integer> edi3alt = new ArrayList(); // creation of an array list to store the alternatives for edi3 by controlling edi1 and edi2
	                for (int foc = 0; foc < edi3opt.size(); foc++) { 
	                	Integer gekor = edi3opt.get(foc);
	                	if (edi1 != gekor && edi2 != gekor && gekor != -1) edi3alt.add(gekor);
	                }
				 	   
	    			// Iterating over the alternatives for the third edge (edi3)
	                for (int telp = 0; telp < edi3alt.size(); telp++) {
	                	int edi3 = edi3alt.get(telp);
	            
	                	NEh = NEs; 
	                	for (int aos = 1; aos < geomit.size(); aos++) {
	                		if(edi3 >= geomit.get(aos-1) && edi3 <= geomit.get(aos)) NEh = geomit.get(aos)-1;
	                	}
	    			
	                	int edip = edi3+1; // ending point id for the third edge delimiting approximate MIC
	                	if (edip >= NEh) edip = edip - NEh;
	                	Coordinate n7 = as[edi3]; // starting point of the third edge
	                	Coordinate n8 = as[edip]; // ending point of the third edge
	            	    LineSegment segd = new LineSegment(n7, n8); // construction of the third edge (d) delimiting MIC
	            
	        	        int tedge = 0; // if the closest point to mic1 or mic2 on the third edge is the ending node of the edge concerned, 
	        	        if (n8.equals(segd.closestPoint(mic1))) tedge = 1; // it is marked with 1 for Apollonius function
	        	  	    if (n8.equals(segd.closestPoint(mic2))) tedge = 1; // it is marked with 1 for Apollonius function
	        	        
	                	// Third approximation for MIC is based on the fact that it should lie at the intersection of the angle bisector of
	        	  	    // the angle formed by the extensions of sega (first edge) and segc (second edge), and the angle bisector of the
	        	  	    // angle formed by the extensions of segc (second edge) and segd (third edge).
	        	        double b44 = segd.angle();
	                    Coordinate teser = segc.lineIntersection(segd); // computation of the intersection of the lines of infinite extent defined by segc and segd.
	                    hacim = 0;
	                    double xp = (n5.x+n7.x)/2; 
	                    double yp = (n5.y+n7.y)/2;
	                    if (teser != null) hacim = Angle.angleBetweenOriented(segc.midPoint(), teser, segd.midPoint())/2;
	                    if (!(Math.abs(hacim) == 0)) { // if the second and third edges (segc and segd) are not parallel to each other
	                    	xp = teser.x;
	                        yp = teser.y;
	                    }
	                	
	                	double an2 = b33 + hacim; // computation of the angle that the vector defined by the angle bisector of segc and segd makes with the X-axis 
	                	if (hacim < 0) an2 = b44 - hacim;
	                	double bdk = Math.tan(an2); // slope of the angle bisector of the lines of segc and segd
	                	double adk = yp - bdk*xp; // intercept of the line represented by the angle bisector of the lines of segc and segd
	                	double x2 = (adk-ack)/(bck-bdk); // the intersection of the angle bisector of sega and segc with the angle bisector of segc and segd
	                	double y2 = ack + bck*x2;
	                	
	                	Coordinate mic3 = new Coordinate(x2, y2); // coordinate of the third approximation for the center of MIC
	                	com.vividsolutions.jts.geom.Point pic3 = gf.createPoint(mic3); // the third approximation for the center of MIC
	                	double rad3 = sega.distance(mic3); // radius of the MIC concerned
	            
	                	if (!pic3.within(geom)) rad3 = 0;
	                	else { // if the third approximation for the MIC overflow the polygon boundaries, the edge closest to mic3 and
	                		bbkng = mic3; // different from the three edges delimiting MIC can also be designated as the third edge. 
	                		bbimp = rad3; // this part also update micCircleRad, micX and micY if the resulting MIC is better than the 
	                		int goler = -1; // existing MIC by re-finding the distance between the closest edge to mic3 and mic3.
	                		for (int lr = 0; lr < doner+1; lr++) { // consideration for polygon with holes ... 
	                			int eki = geomit.get(lr);
	                			LinearRing geomln = geomlin[lr];
	                			Coordinate[] kas = geomln.getCoordinates();
	                			int NEr = kas.length-1;
	                			for (int ip = 0; ip < NEr; ip++) { // finding the closest edge to mic3
	                				double fonlem = CGAlgorithms.distancePointLine(bbkng, kas[ip], kas[ip+1]);
	                				if(fonlem < bbimp) {bbimp = fonlem; goler = ip+eki;}
	                			}
	                		}
		   	   				rad3 = bbimp;
		   	   				if (rad3 >= micCircleRad) {
		   	   					micCircleRad = rad3;
		   	   					micX = mic3.x;
		   	   					micY = mic3.y;
		   	   				}
	            			if (goler != -1 && goler != edi1 && goler != edi2 && goler != edi3) { // if the edge concerned is different
	            				int ekle = 1; // from the three edges delimiting MIC, it can also be considered as an alternative for edi3
	            				for (int fus = 0; fus < edi3alt.size(); fus++) { // provided that it is not already included in edi3alt.
	            					int fusayar = edi3alt.get(fus);
	            					if (goler == fusayar) ekle = 0;
	            				}
	            				if (ekle == 1) edi3alt.add(goler); 
	            			}
	                	}
	            		
	            	    if (segj.distance(mic3) == segj.distancePerpendicular(mic3)){ // if distance between mic3 and segj is equal to perpendicular distance of segj to mic3,
	                    	int jest = 0; // it means that the line passing through mic3 and perpendicular to segj is crossing segj that can be considered as an alternative
	                    	for (int oci = 0; oci < edi2alt.size(); oci++) { // for segc. however, other alternatives for the second edge should be checked in order to  
	                    		double[] juor = edi2alt.get(oci); // control whether or not new alternative is already included in edi2alt.
	                    		int jod2 = (int) juor[0];
	                    		if (jod2 == edis) jest = 1; // starting node of the new alternative edge is edis ...
	                    	}
	                    	if (jest == 0) { // if the new alternative is not included in edi2alt,
	                    		double[] suor = new double[3]; // it is added as an alternative for the second edge into edi2alt.
	                    		suor[0] = edis;
	                    		suor[1] = fuor[1];
	                    		suor[2] = fuor[2];
	     				   		edi2alt.add(suor);
	                    	}
	                    }
	                	
	                	// In this section of the script, spatial configuration of the entities derived from the edges is determined according to four unique combinations
	            	    // In a configuration, there may be (1) three edges (000 - standing for the fact that all entities are edges fedge=0, sedge=0, and tedge=0) that are 
	            	    // already considered above, (2) three points (111 - fedge=1, sedge=1, and tedge=1) for which the best approximation for MIC can already be made by 
	            	    // using the vertices of approximate MA constructed by using Voronoi diagrams, or (3) two edges and one point (001, 010, or 100), or (4) two points 
	            	    // and one edge (011, 101, or 110) that will be addressed by using the specific function created for this purpose: Apollonius.getMIC().
	            
	                	for (int inlp = 0; inlp < edglist.size(); inlp++) { // inner loop created for iteration over the possible configurations of the first and second edges
	                		int [] emor = edglist.get(inlp);
	                		fedge = emor[0];
	                		sedge = emor[1];
	                		
	                		// arrangement of entity configuration into one of the four unique combinations
	              		    ArrayList<Coordinate> taso = new ArrayList<Coordinate>();
	            	        taso.add(n3); taso.add(n4); // 000 into 000, 111 into 111
	            	        taso.add(n5); taso.add(n6); // 001, 010, and 100 into 100
	            	        taso.add(n7); taso.add(n8); // 011, 101, and 110 into 110
	                		int eko = fedge + sedge + tedge;
	                		int ilko = 0; int ikko = 2; int ucko = 4;
	                		if (eko == 1) {
	                			if (fedge == 1) { // 100
	                				ilko = 0; ikko = 2; ucko = 4;
	                			}
	                			else if (sedge == 1) { // 010
	                				ilko = 2; ikko = 0; ucko = 4;
	                			}
	                			else { // 001
	                				ilko = 4; ikko = 0; ucko = 2;
	                			}
	                		}
	                		else if (eko == 2) {
	                			if (fedge == 1) {
	                				ilko = 0;
	                				if (sedge == 1) { // 110
	                					ikko = 2; ucko = 4;
	                				}
	                				else { // 101
	                					ikko = 4; ucko = 2;
	                				}
	                			}
	                			else if (sedge == 1) { // 011
	                				ilko = 2; ikko = 4; ucko = 0;
	                			}
	                		}
				
	                		// for Apollonius.getMIC() there are two cases represented by eko: eko = 1 standing for entity configuration of 001, 010, or 100;
	                		// or eko = 2 standing for entity configuration of 011, 101, or 110. 1 should be arranged in 100 and 2 should be arranged in 110 form.
	                		if (eko == 1 || eko == 2) {
	                			Coordinate api = Apollonius.getMIC(taso.get(ilko+1),taso.get(ikko),taso.get(ikko+1),new LineSegment(taso.get(ucko),taso.get(ucko+1)),eko);
	                			if(api != null) {
	                				pic3 = gf.createPoint(api);
	                				double bic = CGAlgorithms.distancePointLine(api, taso.get(ucko), taso.get(ucko+1));
	                				if (pic3.within(geom)) {
	                					bbkng = api;
	                					bbimp = bic;
	                					for (int lr = 0; lr < doner+1; lr++) { // consideration for polygon with holes ...
	                						LinearRing geomln = geomlin[lr];
	                						Coordinate[] kas = geomln.getCoordinates();
	                						int NEr = kas.length-1;
	                						for (int ip = 0; ip < NEr; ip++) {
	                							double fonlem = CGAlgorithms.distancePointLine(bbkng, kas[ip], kas[ip+1]);
	                							if(fonlem <= bbimp) bbimp = fonlem; 
	                						}
	                					}
	                					if (bbimp >= micCircleRad) {
	                        				micCircleRad = bbimp;
	                        				micX = bbkng.x;
	                        				micY = bbkng.y;
	                					}
	                				}
	                			}
	                		}
		        
	                		// Updating best MIC;         	
	                		if (micCircleRad > micCirclecheck) { 
	                			micCirclecheck = micCircleRad;
	                			micXcheck = micX;
	                			micYcheck = micY;
	                		}
	                		
	                	} // end of loop inlp (inner loop) created for iteration of the possible configurations of the first and second edges
	                } // end of loop telp created for iteration of the alternatives for the third edge 
			 	} // end of loop selp created for iteration of the alternatives for the second edge
			 	   
			 	edi2to1 = (int) edi1alt.get(felp); // as the limited alternatives for the first edge are drawn from the alternatives from the
				if (edi2to1 != -1 && felp == 0) { // second and third edges, they are actually defined in the loop itself and could only be
					ed2i = edi1; // managed in this way ...
				   	edi1 = edi2to1;
				}
				else if (edi2to1 != -1 && felp == 1 && edi2to1 != edi1) {
				   	ed2i = (int) cod[colp][1];
				   	ed3i = edi1;
				   	edi1 = edi2to1;
				}
				else break;
			 } // end of loop felp created for consideration of limited alternatives for the first edge
	  	  } // end of loop colp created for iteration of the core alternatives in the delimitation of MIC
		        
		  micXcheck = micXcheck+precx; // offsetting x coordinate of MIC to its original location
		  micYcheck = micYcheck+precy; // offsetting y coordinate of MIC to its original location
		  
		  com.vividsolutions.jts.geom.Point micnt = gf.createPoint(new Coordinate(micXcheck, micYcheck)); // creating centroid of MIC
		  Geometry ebid = BufferOp.bufferOp(micnt, micCirclecheck, 128); // creating MIC
		  if (!ebid.within(poly)){ // in spite of all controls and checks, the resulting MIC may ignorably overflow or touch the polygon boundary.
		   	   micCirclecheck = micCirclecheck - 0.0000000003; // this may stem from the decimal precision used by the software in the calculations.
		  }
		  
		  double[] mico = new double[3];
		  mico[0] = micXcheck;
		  mico[1] = micYcheck;
		  mico[2] = micCirclecheck;
		  
		  return mico;

	  }
}
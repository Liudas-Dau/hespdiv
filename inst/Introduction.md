---
title: "HespDiv: Concepts and Methodology"
author: "Liudas Daumantas"
date: "2023-06-29"
description: "Vignette provides glossary of key-terms in the HespDiv methodology and the 'hespdiv' package."
output: 
  html_document:
    toc: true
    toc_float: 
      collapsed: false
      smooth_scroll: false
    theme: readable
vignette: >
  %\VignetteIndexEntry{Glossary and Methodology}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## 1. Introduction

Welcome to this introductory tutorial of the "hespdiv" package! This tutorial provides glossary and basic concepts of the methodology implemented in the package. 

The "hespdiv" package provides a specific implementation of the hierarchical spatial data subdivision method called **HespDiv** (not yet published).

Here are the basic characteristics of the HespDiv method:

  1. It is an exploratory spatial data analysis method.
  2. It subdivides data in space in a way that produces subdivisions with an **explicit spatial hierarchical structure**, where each higher-order subdivision is a strict subset of the lower-order subdivision.
  3. The subdivisions produced by the HespDiv method consist of a collection of **split-lines** that divide data in space and   collectively define spatial data clusters (**HespDiv   clusters**) and their spatial boundaries (**HespDiv polygons**).
  4. Each split-line in this output is assigned a **comparison value** that reflects its effectiveness in dividing the data in space.
  5. The HespDiv method is defined primarily by the output it generates and is not tied to a specific implementation algorithm. As a result, there are multiple approaches to conducting HespDiv analysis, and the "hespdiv" package currently provides one of these approaches.
  
Key terms are **bold** in this list, and they will be explained in more detail in the following sections.
  
In addition to the HespDiv implementation algorithm, the "hespdiv"" package provides a range of functions to assess and visualize HespDiv analysis results.

## 2 Overview of The "hespdiv" Package

The "hespdiv" package is designed to enable broad and flexible applications of the HespDiv method, with an emphasis on biogeographical applications. It is not limited to specific data types, research questions, or scientific domains. However, achieving such flexibility requires advanced customization, which may introduce some complexity for users. To address this, the package incorporates several methods specifically tailored for biogeographical studies. These methods use taxa occurrence data to identify hierarchically structured bioregions. While these specialized methods may offer less flexibility in HespDiv applications, their more restrictive nature simplifies the user experience. As the package  evolves, future updates will introduce additional user-friendly options for HespDiv applications. Nonetheless, the overall focus of the package will remain centered on applications in the fields of biogeography, macroecology, and macroevolution.

The "hespdiv" package includes a collection of functions to conduct HespDiv analysis, and assess its results. Additionally, the package is associated with the "HDData" package which provides preloaded data sets of US Miocene land mammal occurrences and US polygon coordinates. These data sets are used for examples throughout the ['HespDiv Walkthrough: Case of US Miocene Mammals'](Walkthrough.Html) vignette that presents a walkthrough case of "hespdiv" package application. 

The main function of the "hespdiv" package is <code>hespdiv()</code>. It implements a somewhat generalized instance of HespDiv algorithm  (indicated by italicized *hespdiv*) that provides several options to modify the optimization algorithm and subdivision criteria. Additionally, four preset methods are integrated in this function that setups the algorithm specifically for bioregionalization purposes.

Other "hespdiv" functions are built to facilitate the interpretation of results. There are three main functions for visualizations: <code>plot_hespdiv()</code>, <code>blok3d()</code>, <code>poly_scheme()</code>. Additionally, there are six functions for *hespdiv* sensitivity analysis: <code>hsa()</code>, <code>hsa_detailed()</code>, <code>change_base()</code> <code>hsa_quant()</code>, <code>plot_hsa_q()</code>, and <code>plot_hsa()</code>. Furthermore, there is one function that allows post-processing of *hespdiv* results and prepares them for analysis (e.g., cluster analysis): <code>cross_comp()</code>. Finally, there is the <code>nulltest()</code> function that allows to test the significance of the obtained split-lines.

## 3. Concepts & Glossary
In this section, you will become familiar with the *hespdiv* algorithm and the specific terminology associated with the HespDiv method and its current implementation.

First, let's clarify the usage of the term 'hespdiv'. The "hespdiv" package includes the <code>hespdiv()</code> function, which implements the *hespdiv* algorithm. It's important to note that while the *hespdiv* algorithm is one approach capable of producing an output characteristic of the HespDiv (refer to Section 1 for the characteristics of the HespDiv method), there is potential for additional implementations of the HespDiv method to be added to the "hespdiv" package in the future.

### 3.1 Glossary
In this subsection, for the purpose of clarity, here is a glossary of key terms used to describe the *hespdiv* algorithm. Refer to this glossary when clarification is needed. The glossary is designed to be read from top to bottom, as terms from above do not require an understanding of terms from below.

Study Area : 
: The study area refers to the geographic region or space defined by a polygon. By default this polygon is obtained by calculating a convex hull of observation locations (<code>**'xy.dat'**</code>). However, it can also be provided directly in the input (<code>**'study.pol'**</code>). When the algorithm description mentions "*study area*", it refers to the original study area polygon (or region itself) and not its subsets obtained through subdivisions.

Study Data : 
: Similar to the study area, study area refers to the original input data. It contains information about observations (<code>**'data'**</code>) and their locations (<code>**'xy.dat'**</code>).

Subdivision : 
: Subdivision can either refer to the process of subdividing the study area and data or the result of such a process. The interpretation depends on the context.

HespDiv polygon : 
: The first HespDiv polygon is the polygon representing the study area. Further HespDiv polygons are created by subdividing this polygon.

A HespDiv cluster : 
: A HespDiv cluster is a spatial cluster of study data points obtained by filtering the study data with a corresponding HespDiv polygon. The first HespDiv cluster encompasses all study data.

Split-lines :  
: Split-lines are lines that divide HespDiv polygons and clusters into two parts. They can be either straight/linear or nonlinear.

Rank (or order) of a HespDiv polygon or cluster : 
: The Rank of HespDiv polygon or cluster refers to its location in the spatial hierarchy tree. The HespDiv polygon and cluster of 1st rank/order encompass all study area and data, respectively. Subsequent subdivisions produce higher order HespDiv polygons and clusters that are strict subsets of lower order HespDiv polygons and clusters, respectively. For instance, if a HespDiv polygon or cluster is of rank five then it means that it took four subdivisions to produce it.

Split-points :  
: Split-points are evenly distributed points on the perimeter of a HespDiv polygon. Connecting them creates straight/linear split-lines.

Straight/linear split-lines:
: A straight/linear split-line is created by connecting two split-points with a straight line.

Subdivision Method :
: The subdivision method determines how the quality of split-lines is measured and optimized. 

  : It consists of three steps: 

1) Obtaining HespDiv polygon objects by transforming or generalizing each HespDiv cluster from opposite sides of a split-line  (*using* <code>**'generalize.f()'**</code>).
2) Comparing these objects to obtain a comparison value for a split-line (*using* <code>**'compare.f()'**</code>).
3) Selecting split-lines that produce higher or lower comparison values (*using* <code>**'maximize'**</code>).

HespDiv polygon object : 
: A HespDiv polygon object is obtained after applying the first step of the subdivision method to a HespDiv cluster. It represents a transformed or generalized subset of study data, and its interpretation depends on the used subdivision method (e.g., a list of unique taxa names or a table of taxa composition).

Comparison Value : 
: A comparison value is used to describe the strength or weakness of a subdivision produced by a split-line. It is obtained in the second step of the subdivision method by comparing HespDiv polygon objects obtained from opposite sides of a split-line. The interpretation of the comparison value is defined by the subdivision method used (e.g., Morisita-Horn similarity index).

Split-line Performance/Quality:
: The split-line performance or quality terms in this text do not depend on the metric of comparison values. Thus, if a low comparison value indicates good subdivision, the performance or quality of the split-line producing this subdivision should be interpreted as an inverse of its comparison value.

The Best Split-line :  
: The best split-line, whether linear or nonlinear, represents the split-line with the highest performance.

Subdivision criteria:
: Subdivision criteria are conditions that split-lines and the resulting HespDiv polygons and clusters must meet to be considered for subdivision.

Absolute criteria:
: Absolute criteria are subdivision criteria whose values should be interpreted with respect to the study area and study data (the first-order HespDiv polygon and cluster).

Relative criteria:
: Relative criteria are subdivision criteria whose values should be interpreted with respect to the HespDiv polygon and cluster currently being subdivided. Therefore, relative criteria are always provided as proportions.

Location number criterion:
: The location number criterion specifies the minimum number of observation locations that must be in each HespDiv polygon produced by a split-line for it to be considered for subdivision. This criterion has absolute and relative versions.

Area size criterion:
: The area size criterion specifies the minimum area of each HespDiv polygon produced by a split-line for it to be considered for subdivision. This criterion has absolute and relative versions.

Sample size criterion:
: The sample size criterion specifies the minimum number of observations that must be in each HespDiv polygon produced by a split-line for it to be considered for subdivision. This criterion has absolute and relative versions.

Quality criterion :
: A quality criterion is a minimum performance that a split-line must achieve to be selected for subdivision.

Testing a split-line:
: The process of testing a split-line, whether linear or nonlinear, involves one or three steps: 

  1) Checking if the two HespDiv polygons and clusters obtained with a split-line division meet the relative and absolute criteria set for area size, location number, and sample size. If not, testing halts.
  2) Applying a subdivision method to obtain a comparison value for a split-line.
  3) Comparing the obtained performance value to a performance threshold value to determine if it surpasses it. Initially, the performance threshold is the provided quality criterion from the input, but it later becomes the comparison value of the best split-line that has surpassed the quality criterion. If a split-line surpasses the performance threshold, it is selected for subdivision unless a better split-line is found in subsequent tests.

Nonlinear split-lines :  
: Nonlinear split-lines deviate from the best straight split-line and are generated using smooth mathematical functions called splines to divide HespDiv polygons and clusters.

A Network (or Net) of (Spline) Knots / A (Spline) Knot Network : 
: A net of spline knots refers to a set of control points allocated around the best straight split-line. These knots serve as input for generating nonlinear split-lines using spline interpolation. The knots in this network are distributed systematically, with knots placed at regular intervals on straight lines orthogonal to the best straight split-line. These lines are parallel to each other and evenly spaced along the best straight split-line.

Column of knots (or knot column):
: A knot column is a set of knots in the knot network aligned on a straight line orthogonal to the best straight split-line. Each knot column contains the same number of evenly distributed knots, but the spacing between knots may vary based on the HespDiv polygon's shape.

Knot combination:
: In essence, different knot combinations can be understood as different nonlinear split-lines, as each knot combination, when interpolated using spline interpolation, produces a distinct nonlinear split-line. Each knot combination contains one knot from each knot column, resulting in as many knots as there are knot columns. Therefore, the number of knot columns determines the possible number of wiggles in the generated split-lines.

Testing a knot combination:
: The process of testing a knot combination involves using spline interpolation to generate a nonlinear split-line from the knot combination and testing the obtained split-line.

Testing a knot:
: The process of testing a knot involves obtaining a knot combination by replacing a knot in the best knot combination with a tested knot and testing the obtained knot combination. The replaced knot is always from the same knot column as the knot being tested. In a sense, testing a knot means the same as testing a nonlinear split-line.

The Best Combination of (Spline) Knots / The Best (Spline) Knot Combination :  
: The best combination of (spline) knots is a specific combination of knots that, when interpolated using spline interpolation, produces the best nonlinear split-line among all tested combinations. Each knot in this combination is from a different knot column. Initially, it is assumed that the best combination of knots lies directly on the best straight split-line, where lines corresponding to knot columns intersect it. When interpolated, this particular set of knots produces a split-line that is identical to the best straight split-line.

Iteration through the knot network:
: The hespdiv algorithm repeatedly iterates through the knot network to find the best knot combination. It does so by iterating through knot columns forth and back along the best straight split-line, always testing different combinations of knots. When iterating through a particular knot column, the algorithm iterates through its knots, testing each one of them. Iteration stops when no new knot combinations can be obtained. Alternatively, the default behavior can be changed by setting a limit on the maximum number of iterations through the knot columns (*using* <code>**'c.max.iter.no'**</code>).

Knot performance:  
: Knot performance refers to the performance of a nonlinear split-line obtained during the process of testing a knot.

The best knot:
: The best knot is the knot with the highest performance.

Knot selection:
: Knot selection is the process of updating the best knot combination with the best knot. If a knot test does not yield the highest performance, the best knot combination remains the same as before the test. Otherwise, the default behavior of the algorithm would be to finish the ongoing iteration through knot columns and then replace the current best combination of knots with the knot combination that includes the best knot. Alternatively, this replacement can be performed as soon as the best knot is found (*using* <code>**'c.fast.optim'**</code>). This way, computation would be a bit faster at a small cost to performance.

Interpolated (Spline) Knots :  
: Interpolated (spline) knots are additional knots generated and tested when the algorithm finishes iterating through the knots of a single knot column. The locations of these knots are determined by interpolating knot performance as a function of location in a knot column. An interpolated knot is placed at a location in the knot column where the maximum interpolated knot performance is estimated. Interpolated knots are tested in the same manner as regular knots. If selected, they are included in the best knot combination but not added to the knot network.

## 4. The Methodology

In essence, the *hespdiv* algorithm systematically generates and tests a number of linear and nonlinear split-lines to divide the study data and study area into subsets in order to find the optimal subdivision. This process is repeated recursively with the obtained subsets until no further subdivision can be established due to unmet subdivision criteria.

A more detailed description of the default *hespdiv* behavior would be as follows: 

The algorithm creates a convex hull of coordinates of provided observation locations. This convex hull serves as the study area and the first HespDiv polygon. Then the algorithm starts the recursive iterations. Each iteration it starts by generating straight split-lines (also called linear split-lines). It does so by placing a predetermined number of split-points along the HespDiv polygon and by connecting them with straight lines. Each linear split-line is then tested in a manner that depends on the chosen subdivision method, aiming to identify the one that produces the most optimal subdivision. The best linear split-line is then used as a basis to generate curves (nonlinear split-lines). These curves are created using smooth mathematical functions called splines. A network of spline knots (control points for curve shape) is allocated around the best straight split-line. At first, it is assumed that the best combination of knots is composed of the knots that are placed directly on the best straight split-line. This particular knot combination produces a split-line that is equal to the best straight split-line. The algorithm then iterates through the knot network, at a time making one knot change in the best knot combination, each time producing a curve that diverges from the best straight split-line at the location of the added knot. These curves are tested in the same manner as straight split-lines. When the algorithm finishes iterating through all knot network, it makes one change to the assumed best knot combination, adding the best knot. Then it iterates again and again through the knot network until knot combination adjustments no more results in curves of better performance. Finally, if the best nonlinear split-line supersedes the best straight split-line by a desired margin, it is selected to perform subdivision. Otherwise, subdivision is performed with the best straight split-line. The subdivision ends the current recursive iteration and the obtained two subsets of data and study area are further subdivided in subsequent recursive iterations.

## 5. Algorithm & Arguments

The number of arguments in <code>hespdiv()</code> function may seem overwhelming at first. However, for a basic usage, only a small subset of these arguments is necessary, as default values are selected to suit many cases. Other arguments are dedicated to advanced usage only. 

The main groups of arguments, based on their influence or purpose, are as follows:

  - Study data: *'data'* and *'xy.dat'*.
  - Study area: *'study.pol'* or *'use.chull'*.
  - Subdivision method:
    - Basic usage: *'method'*.
    - Advanced usage: *'compare.f'*, *'generalize.f'*, *'maximize'*.
  - Subdivision criteria:
    - Sample size criteria: *'N.crit'*, *'N.rel.crit'*.
    - Location number criteria: *'N.loc.crit'*, *'N.rel.loc.crit'*.
    - Area criteria: *'S.crit'*, *'S.rel.crit'*.
    - Quality criteria: *'Q.crit'* and *'c.crit.improv'* (also *'c.Q.crit'* in advanced cases).
  - Fit to data and computation speed:
    - Straight split-lines: *'n.split.pts'* (also influenced by *'same.n.split'*).
    - Nonlinear split-lines: *'c.X.knots'* and *'c.X.knots'* (also influenced by *'c.max.iter.no'* and *'c.fast.optim'*).
  - Algorithm customization (advanced usage):
    - Straight split-lines: *'same.n.split'*.
    - Type of split-lines: *'c.splits'*.
    - Nonlinear split-lines: *'c.max.iter.no'* and *'c.fast.optim'*.
  - Miscellaneous: *'tracing'*, *'pnts.col'*, *'display'*, *'c.corr.term'*, *'pacific.region'*.
  

Now let's take a closer look at how the *hespdiv* algorithm uses these arguments to obtain hierarchical spatial data subdivisions:

  1) **Input Pre-processing:**
  - Study Data: Study data is defined by two arguments ('data' and 'xy.dat') to increase the flexibility of *hespdiv*. The rows of 'xy.dat' data frame provides coordinates for the observations (rows or elements) in the 'data' object. The structure of the 'data' depends on the subdivision method. The currently integrated methods require 'data' to be a numeric or character vector of taxa indexes or names, respectively. 
   - Study Area: By default, algorithm obtains a convex hull from the 'xy.dat' data frame and uses it as a study area and the first HespDiv polygon. However, you can provide a data frame of polygon vertex coordinates using the 'study.pol' argument and set the 'use.chull' argument to 'FALSE' to use this polygon as a study area and the first HespDiv polygon.  In specific cases where the study area straddles the 180th meridian, you should set the 'pacific.region' to 'TRUE'. This will transform 'xy.dat' and 'study.pol' (if provided) in a way that enables the analysis to be performed.
   - Subdivision Method: Subdivision method controls how the quality of split-lines is measured and optimized. The method is defined by a combination of three arguments: compare.f, generalize.f and maximize. These arguments are defined internally if you select a preset subdivision method. There are four preset options: 'pielou', 'morisita', 'sorensen', 'horn.morisita'. All of them are designed for hierarchical bioregionalization. For more information about these methods, refer to the documentation of the <code>hespdiv()</code> function.

2) **Split-point Placement:** The algorithm begins by evenly placing a predetermined number of split-points along the perimeter of a HespDiv polygon. The number of points placed is determined by <code>n.split.pts + 1</code>. The 'same.n.split' argument also impacts the quantity of these points. When set to 'TRUE', the number of points placed remains the same regardless of the polygon size. When set to 'FALSE', the spacing between points remains consistent, leading to fewer points for smaller polygons created by higher order subdivisions. By default 16 split-points are
placed on the perimeter, which results in a distance between split-points equal
to 1/16th of the polygon circumference. 

3) **Straight Split-lines:** Straight split-lines are generated by connecting the split-points. The total number of straight split-lines generated is determined by the formula <code>sum(1:n.split.pts)</code>. However, this formula applies only in the first iteration when 'same.n.split' is set to 'FALSE'. It's important to note that not all generated split-lines are evaluated. Split-lines that cross polygon boundaries or fail to meet the sample size, area, or location number criteria are not evaluated.

4) **Subdivisions:** Each straight split-line that does not cross polygon boundaries spatially divides the study data and HespDiv polygon into two subsets.

5) **Criteria:** Both subsets are then checked to determine if they meet sample size, area and location number absolute and relative subdivision criteria.

6) **Comparison Values**: Subsets that meet the criteria are compared using generalize.f and compare.f functions to obtain a comparison value. These functions are set internally, if an integrated subdivision method is used. First, each subset is passed into the generalize.f function to obtain a generalization value (e.g., Pielou evenness index in 'pielou' method) for it. Then these values are used in the compare.f function to compare the subsets, producing a comparison value. This way each split-line that passed all subdivision criteria is assigned a comparison value.

7) **Straight Split-line Selection**: The best performing straight split-line is determined based on whether the Boolean value of maximize argument. This argument is also defined internally, if integrated subdivision method is used. If maximize is TRUE, the best split-line is the one with the highest comparison value; if maximize is FALSE, the best split-line has the lowest comparison value.

8) **Nonlinear Split-lines:** The default behavior of the algorithm is to attempt to increase the performance of the best straight split-line by bending it, giving it a nonlinear shape. This behavior can be changed using the 'c.splits' argument.

If the best straight split-line meets the quality criteria specified by the 'c.Q.crit' argument, it serves as the basis for generating variously shaped curves (nonlinear split-lines). These curves are produced using splines, which are mathematical functions that can create smooth and flexible curves.

To generate the curves, a number of knots (control points) for the splines are distributed evenly along a set of lines orthogonal to the best straight split-line. These orthogonal lines are also evenly distributed along the straight split-line itself. The number of knots and lines used can be adjusted through the arguments 'c.Y.knots' and 'c.X.knots', respectively. 

The algorithm iterates through this network of knots, considering different combinations of knots, to produce curves. By varying the selection and arrangement of knots, different shapes of curves are generated. Adding more lines orthogonal to the straight split-line with the 'c.X.knots' argument allows more wiggles in the produced curves, whereas adding more knots in each line with the 'c.Y.knots' argument allows each wiggle to achieve more shapes.

9) **Nonlinear Split-line Evaluation and Selection**: Curves are then processed in the same manner as straight split-lines in steps 4 to 7.

10) **Final Split-line Selection and Establishment of Subdivision:** If the best curve outperforms the best straight split-line by a margin of c.crit.improv, the best split-line becomes nonlinear. Otherwise, it remains straight. The split-line must also satisfy the criteria established by the Q.crit argument to be used for subdivision.

11) **Recursive Iteration:** The process described from step 2 to 10 is iteratively applied, resulting in a collection of selected split-lines with their performance values. These split-lines hierarchically subdivide space and data, forming polygons of various shapes.

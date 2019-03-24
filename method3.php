<?php


//significance threshold for adjusted p-values in the 4 differential expression analyses
$signif = 0.01;

//River Plate score threshold parameter (R)
$threshold = 4.0;

//whether to turn zeroes into ones and vice-versa, for HyperTraPS
$invert = 0;

//whether to intercalate extra lines of "final states", for HyperTraPS
$addextraline = 0;

//whether to include the sample names in the generated output
$outputsamplenames = 1;

//whether to include the parameter values in the generated output
$outputspecs = 1;

//prefix for the file containing the generated outout
$prefix = "matrix";
//$prefix = "hypertraps";


//------------------------------------------------------------------------------------------------------------
$num_genes = 54675; //number of genes (probes) in the dataset
$sampletypes = "00001111111111112222222222222222233333333333333333333444444444444444";  //4, 12, 17, 20 and 15 in each class
//------------------------------------------------------------------------------------------------------------


//table of what genes are differentially expressed in each stage vs control
$de = array();

//see what stage(s) each gene belongs to
for ($stage = 1 ; $stage<=4 ; $stage++) {
	$dea = explode("\n",file_get_contents($stage."vsC.csv"));
	for ($gene = 1 ; $gene <= $num_genes ; $gene++) {
		$line = explode(",",$dea[$gene]);
		if ($line[5] < $signif) $de[$stage][rq($line[0])] = 1;
	}
}

//open matrix
$ex = explode("\n",file_get_contents("ex.csv"));

//get sample names from first line of "ex"
$samplenames = array();
$line = explode(",",$ex[0]);
for ($sample = 1 ; $sample <= 68 ; $sample++) $samplenames[$sample] = rq($line[$sample]);

//initialise arrays for the results
$squares=array();
$belongs=array();
$squareexps=array();

//loop through the genes
for ($gene = 1 ; $gene <= $num_genes ; $gene++) {
	
	//a line from the ex file
	$line = explode(",",$ex[$gene]);
	
	//extract the gene id from that line
	$geneid = rq($line[0]);
	
	//decide what square this gene belongs to
	$square = "";
	for ($stage = 1 ; $stage<=4 ; $stage++) $square .= isset($de[$stage][$geneid]) ? "Y" : "N";
	if ($square == "NNNN") continue;
	
	//add this gene to the correct square
	$squares[$square] []= $geneid;
	
	//record this gene as belonging to a certain square
	$belongs[$geneid] = $square;
	
	//remember the f
	for ($sample = 1 ; $sample <=68 ; $sample++) {
			$stage = $sampletypes[$sample-1];
			if (!isset($squareexps[$square][$geneid][$stage])) $squareexps[$square][$geneid][$stage]=0;
			$squareexps[$square][$geneid][$stage] += $line[$sample];
	}
	
}

//divide the sums of squares by the number of items in each subset to obtain averages
foreach ($squareexps as $square => $genes) {
	foreach ($genes as $gene => $values) {
		//4 12 17 20 15
		$squareexps[$square][$gene][0] /= 4;
		$squareexps[$square][$gene][1] /= 12;
		$squareexps[$square][$gene][2] /= 17;
		$squareexps[$square][$gene][3] /= 20;
		$squareexps[$square][$gene][4] /= 15;
	}
}
	

//the 16 subsets
$features = array("NNNN","NNNY","NNYN","NNYY","NYNN","NYNY","NYYN","NYYY","YNNN","YNNY","YNYN","YNYY","YYNN","YYNY","YYYN","YYYY");

//the 16 subsets minus the discarded "NNNN" subset
$rfeatures = array(      "NNNY","NNYN","NNYY","NYNN","NYNY","NYYN","NYYY","YNNN","YNNY","YNYN","YNYY","YYNN","YYNY","YYYN","YYYY");

//load, modify and output the table of groups
$magic = file_get_contents("magic");
foreach ($features as $feature) {
	$square = "";
	$square .= $feature . "<br>";
	$square .= count($squares[$feature]);
	$magic = str_replace("c$feature",$square,$magic);
}
echo $magic;

//initialise scores matrix
$score = array();
for ($sample = 1 ; $sample <= 68 ; $sample++) {
	foreach ($rfeatures as $feature) {
		$score[$sample][$feature] = 0;
	}
}

//initialise output
$hypertraps = "";

//output the values of the parameters used
$specs = "Method = 3 ; Threshold = $threshold";
echo "<h2>$specs</h2>";
if ($outputspecs) $hypertraps .= "$specs\n";


//array to store the expression levels for a sample
$sexp = array();


//output table header
echo "<table border=1>";
echo "<tr><th>Feature</th>";
foreach ($rfeatures as $feature) echo "<th>$feature</th>";
echo "</tr><tr><th># genes</th>";
foreach ($rfeatures as $feature) echo "<th>".count($squareexps[$feature])."</th>";


//determines ones and zeroes for each sample
for ($sample = 1 ; $sample <=68 ; $sample++) {
	
	//loop through the genes
	for ($gene = 1 ; $gene <= $num_genes ; $gene++) {
		$line = explode(",",$ex[$gene]);
		$geneid = rq($line[0]);
		$sexp[$geneid] = $line[$sample];
	}
	
	$addtoht = array(); //items to add to the file for hypertraps
	$scores = array();
	
	//loop through the subsets (groups)
	foreach ($rfeatures as $feature) {
	
		foreach ($squareexps[$feature] as $geneid => $avgs) {
			$avg = $avgs[0];
			$score[$sample][$feature] +=  ($sexp[$geneid]-$avg)*($sexp[$geneid]-$avg);
		}
		$score[$sample][$feature] /= count($squareexps[$feature]); //score for this sample for this group of genes
		
		$bit = ($score[$sample][$feature] > $threshold) ? 1 : 0;  //determine if it's a 1 or a 0
			
		if ($invert) $bit = 1-$bit; //invert if required for hypertraps
		$addtoht []= $bit;
		$scores []= sprintf("%.2f",$score[$sample][$feature]);
	}
	
	if ($outputsamplenames) $hypertraps .= $samplenames[$sample]." "; //output sample names if required
	
	$hypertraps .= implode(" ",$addtoht);
	$hypertraps .= "\n";
	if ($addextraline) { //add extra line of "final states" if required for hypertraps
		$hypertraps .= "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
		$hypertraps .= "\n";
	}
	echo "<tr><td style='background-color:".(($sampletypes[$sample-1]%2)?"#DDDDDD;":"white")."'>$samplenames[$sample] (<b>".$sampletypes[$sample-1]."</b>)</td>";
	foreach ($scores as $s) echo "<td style='text-align:center;color:white;background-color:".(($s>$threshold)?"red":"green").";'><b>$s</b></td>";
	echo "</tr>";
	
}

echo "</table>";

echo "<hr>";

echo "<pre>";
echo $hypertraps;
echo "</pre>";

file_put_contents($prefix.$signif."-".$threshold.".txt",$hypertraps); //save file for hypertraps



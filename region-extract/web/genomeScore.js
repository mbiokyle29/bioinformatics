
var fastaHeader;
var sequence;
var scaler;

d3.text("test.fasta", function(sequenceData) {
	var splitData = sequenceData.split("\n");
	
	fastaHeader = splitData[0];
	sequence = splitData[1].split("");
	scaler = d3.scale.linear()
		.domain([0, sequence.length])
		.range([0, sequence.length]);

	var container = d3.select(".d3-box")
	.append("svg")
    	.attr("width", 1000) // scaler(sequence.length))
    	.attr("height", 1000);

    var yPos = 100;
    var xPos = 10;

    container.selectAll("text")
    	.data(sequence)
    	.enter()
    	.append("text")
    	.attr('y', yPos)
    	.attr('x', function(d) { xPos += 15; return xPos; })
    	.text(function(d) {  return d; })
    		.attr("font-family", "sans-serif")
			.attr("font-size", "20px")
			.attr("fill", "red");

	//var bases = container.selectAll("text")
	//	.data(sequence)
	//	.enter()
	//	.append("text")
	//	.text(function(d) { return d; });

});
var test = ['A', 'T', 'C', 'G'];



var posScore = d3.text("kozak.plus.out");
var minScore = d3.text("kozak.minus.out");




# Define a Python subroutine to colour atoms by B-factor, using predefined
# intervals


def colour_consurf(selection="all"):

	# These are constants
	min = 0.0
	max = 9.0
	n_colours = 10
	colours = [ 
	[1, 1, 0.58823529] ,
	[0.062745098,0.78431373,0.81960784], 
	[0.54901961,1,1], 
	[0.84313725,1,1], 
	[0.91764706,1,1], 
	[1,1,1], 
	[0.98823529,0.92941176,0.95686275], 
	[0.98039216,0.78823529,0.87058824], 
	[0.94117647,0.49019608,0.67058824], 
	[0.62745098,0.14509804,0.37647059] ]
	
	bin_size = ((max - min) + 1) / n_colours

	# Loop through colour intervals
	for i in range(n_colours):

		lower = min + i * bin_size  
		upper = lower + bin_size    
		colour = colours[i]

		# Print out B-factor limits and the colour for this group
		print lower, " - ", upper, " = ", colour

		# Define a unique name for the atoms which fall into this group
		group = selection + "_group_" + str(i+1)
		
		# Compose a selection command which will select all atoms which are 
		#	a) in the original selection, AND
		#	b) have B factor in range lower <= b < upper
		sel_string = selection + " & ! b < " + str(lower)
		
		if(i < n_colours - 1):
			sel_string += " & b < " + str(upper)
		else:
			sel_string += " & ! b > " + str(upper)
		
		# Select the atoms
		cmd.select(group, sel_string) 

		# Create a new colour
		colour_name = "colour_" + str(i+1)
		cmd.set_color(colour_name, colour)

		# Colour them
		cmd.color(colour_name, group)


	# Create new colour for insufficient sequences
	insuf_colour = [1, 1, 0.58823529]
	cmd.set_color("insufficient_colour", insuf_colour)

	# Colour atoms with B-factor of 10 using the new colour
	cmd.select("insufficient", selection + " & b = 10")
	cmd.color("insufficient_colour", "insufficient")




# This is required to make command available in PyMOL 
cmd.extend("colour_consurf", colour_consurf)

colour_consurf()

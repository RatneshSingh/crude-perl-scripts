my$num_colors=$ARGV[0];
$num_colors=$num_colors?$num_colors:10;
open OUT,">Color.html";
print "\n";
print OUT "\<\!DOCTYPE html\>\n\<html\>\n\<body\>\n";
foreach my$color(sort(colourset($num_colors))){
        my$rgb=$color;
        $rgb=~s/^\s+//g;
        $rgb=~s/\s+$//g;

        $rgb=~s/\s+/,/g;
    #print OUT "\n\<p style\=\"color:rgb\($rgb\)\"\>$color<\/p\>\n";
        print "$color\n";
   print OUT "\<div style\=\"position:relative\;left:50px\;width:300px\;height:50px\;text-align: center\;background-color:rgb\($rgb\)\"\>$color\<\/div\>\n";
}

print OUT "\<\/body\>\<\/html\>";
close OUT;
print "\n\n\nOne Html file named Color.html has been saved in current folder. Please use it to visualize the colors\n\n";


#####sub routines
#=item colourset($num_colours, $method)
#Function to grab a set of well spaced colours from an HSV colour wheel.
#    $num_colours - The number of colour values to produce, must be greater than 0 but no bigger than 360
#    $method - The method for selecting colours over HSV colour space, either 'equal_spacing' or for around 10 colours 'chroma_bisection' is better.
#Returns an array of RGB values of the form ([R,G,B]) and undef on $num_colours out of bounds
#=cut
sub colourset {
        my ($num_colours,$method) = @_;
        if ($num_colours <= 0 or $num_colours > 360) {
                warn "Number of colours requested out of bounds.";
                return undef;
        }
        $method = 'chroma_bisection' unless $method;

        ##Internal sub to randomly shuffle an array
        #sub fisher_yates_shuffle {
        #        my ($array) = @_;
        #        my $current;
        #        for ($current = @$array; --$current; ) {
        #                my $selected = int rand ($current+1);
        #                next if $current == $selected;
        #                #Reverse the slice between current position and the randomly selected
        #                @$array[$current,$selected] = @$array[$selected,$current];
        #        }
        #        return $array;
        #}

        #Colours to return
        my %colours;

        #Default Hue Saturation and Value, saturation of 0.65 gives a more pastel feel!
        my ($Hue, $Saturation, $Value) = (0.0,0.65,1.0);

        #The interval to space colours around the wheel if equal
        my $hsv_interval = 360 / $num_colours;

        #Array of degrees for reuse to create ranged arrays with a given interval
        my @degrees = 1..360;

        #Iteratively bisect each chroma segment so that the first 6 colours are well spaced perceptually.
        #However after 12 colours we will have increasing pairs that are more confused as
        #they are increasingly close to each other compared to the rest of the colours!
        #To get around this problem of selecting closely around a single bisection, we jump around the
        #chroma randomly sampling.
        if ($method eq 'chroma_bisection') {
                #The current cycle of chroma bisection
                my $hsv_cycle = 1;
                #Number of colours selected by bisecting chroma so far
                my $colours_selected = 0;

                #While we still have colours to select
                while ($colours_selected != $num_colours) {
                        #Work out the size of interval to use this cycle around the wheel
                        $hsv_interval = 60 / $hsv_cycle;
                        #Get all the Hues for this cycle that haven't already been examined and are on the line of bisection
                        my @Hues = grep { (not $_ % $hsv_interval) && (not exists $colours{$_%360}) } @degrees;
                        #Shuffle so that we don't take from around the same chroma all the time, only perceptually worthwhile after 12th colour
                        fisher_yates_shuffle(\@Hues) if $hsv_cycle > 2;

                        #While we still have hues to select from in this cycle
                        while (@Hues) {
                                #Finish if we have enough colours
                                last if $colours_selected == $num_colours;
                                #Consume a Hue from this cycle
                                $Hue = shift @Hues;
                                #360 should be 0 for red
                                $Hue %= 360;
                                #Store this Hue and mark selection
                                $colours{$Hue} = hsv2rgb($Hue,$Saturation,$Value) ;
                                #print "\n*****RGB value of color: $colours{$Hue}";
                                $colours_selected++;
                        }
                        $hsv_cycle++;
                }
        }

        #Just space colours even distances apart over the HSV colour wheel.
        #You have slightly odd/garish colours coming out, but you dont get uneven perceptual distance
        #between pairs of colours. This scales far better despite the horrible colours.
        elsif ($method eq 'equal_spacing') {
                foreach $Hue (1..$num_colours) {
                        $Hue = ($Hue * $hsv_interval) % 360;
                        $colours{$Hue} = hsv2rgb($Hue,$Saturation,$Value) ;
                }
        }

        #Otherwise return nothing and warn the programmer
        else {
                warn "Colourset method not known, use either 'equal_spacing' or for fewer colours 'chroma_bisection'";
                return undef;
        }

        #Shuffle final colours so that even if we do use chroma_bisection closer colours will hopefully not be sequential
        @_ = values %colours;
        fisher_yates_shuffle(\@_);
        return @_;
}

######## HSV to RGB in Perl
#
######## And for completeness here is the Perl function to convert from HSV to RGB colour space, essentially this is a mapping from a cylinder to a cube with a rotation. I can't claim to have worked all this out, I just ported a function I found online.
#
#=item hsv2rgb($Hue, $Saturation, $Value)
#Function to convert HSV colour space values to RGB colour space.
#Returns RGB value as [R,G,B]
#=cut
sub hsv2rgb {
        my ($Hue,$Saturation,$Value) = @_;
        my ($Red,$Green,$Blue) = (0,0,0);

        #Check the input and warn if it's a bit wrong
        warn "Invalid Hue component of HSV colour passed, with value: $Hue." unless ($Hue >= 0.0 and $Hue <= 360.0);
        warn "Invalid Saturation component of HSV colour passed, with value: $Saturation." unless($Saturation >= 0.0 and $Saturation <= 1.0);
        warn "Invalid Value component of HSV colour passed, with value: $Value." unless ($Value >= 0.0 and $Value <= 1.0);

        #If colour has no saturation return greyscale RGB
        if ($Saturation == 0) {
                $Red = $Green = $Blue = $Value;
                return [$Red, $Green, $Blue];
        }

        #Partition the Hue into the 5 different colour chroma and then map each of these to RGB based on the colour theory
        $Hue /= 60.0;
        use POSIX;
        my $Chroma = floor($Hue) % 6;
        my $H_d = $Hue - $Chroma;

        #RGB cube components
        my ($I,$J,$K) = ( $Value * ( 1 - $Saturation ),
                                   $Value * ( 1 - $Saturation * $H_d ),
                                   $Value * ( 1 - $Saturation * ( 1 - $H_d ) )
                                    );

        #Map components to RGB values per chroma
        if ($Chroma == 0) { ($Red,$Green,$Blue) = ($Value,$K,$I); }
        elsif ($Chroma == 1) { ($Red,$Green,$Blue) = ($J,$Value,$I); }
        elsif ($Chroma == 2) { ($Red,$Green,$Blue) = ($I,$Value,$K); }
        elsif ($Chroma == 3) { ($Red,$Green,$Blue) = ($I,$J,$Value); }
        elsif ($Chroma == 4) { ($Red,$Green,$Blue) = ($K,$I,$Value); }
        else{ ($Red,$Green,$Blue) = ($Value,$I,$J); }
        #print "\nReturning Red:$Red,Green:$Green,Blue:$Blue";
        #Return the RGB value in the integer range [0,255] rather than real [0,1]
        #return [floor($Red * 255),floor($Green * 255),floor($Blue * 255)];
        return join(" ",floor($Red * 255), floor($Green * 255), floor($Blue * 255));
}


sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = scalar@$array - 1;$i>=0; --$i ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
    #print "\n******colors came to shuffle $array";
}
#!/usr/bin/awk

BEGIN{
    print "Chromosome\tPosition\tReference\tAlternate\tQuality\tDepth\tAlt_Freq\tAnnotation";
}
$5 ~ /[ATCG],*[ATCG]*/ {
    if($1 !~ /#.*/){
        a=split($8, Dat, ";");
        #for (key in Dat){ print Dat[key] }
        split(Dat[4], Ann, "=");
        split(Dat[1], AF, "=")
        split(Dat[3], DPt, "=")
        if(AF[2] > f && DPt[2] > d){
            print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" DPt[2] "\t" AF[2] "\t" Ann[2];
        }
    }
}
END{}
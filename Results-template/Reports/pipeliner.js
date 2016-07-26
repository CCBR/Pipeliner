function fastqc() {
    var x = document.fastqcform.fastqcselect.value;
    var y="../QC/"+x+".R1_fastqc.html";
    var z="../QC/"+x+".R2_fastqc.html";
    document.fastqcform.fastqcr1html.data=y;
    document.fastqcform.fastqcr2html.data=z;
//    alert("name= "+y);
}

function fastqctrimmed() {
    var x = document.fastqcform.fastqcselecttrimmed.value;
    var y="../QC/"+x+".R1.trimmed_fastqc.html";
    var z="../QC/"+x+".R2.trimmed_fastqc.html";
    document.fastqcform.fastqcr1trimmedhtml.data=y;
    document.fastqcform.fastqcr2trimmedhtml.data=z;
//    alert("name= "+y);
}

function qualimap() {
    var x = document.fastqcform.qualimapselect.value;
    var y="../QC/"+x+".qualimapReport/qualimapReport.html";
    document.fastqcform.qualimaphtml.data=y;

}

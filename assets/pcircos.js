alert('If you see this alert, then your custom JavaScript script has run!')
$(document).ready(function() {
    alert("document is ready!")
    
    $( ".sidebar" ).resizable({
        handles: "e, w"
    });
    $( ".graph" ).resizable({
        handles: "e, w"
    });
    $( "summary" ).click(function () {
        $( this ).toggleClass("toggle")
    })
}
)

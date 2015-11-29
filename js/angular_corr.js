function plot(){
    var graph = document.getElementById("graph_div");
    var width = parseInt(graph.style.width, 10);
    var x1 = -1;
    var x2 = 1;
    var a2 = parseFloat(document.getElementById("a2").value);
    var a4 = parseFloat(document.getElementById("a4").value);
    var xs = 1.0 * (x2 - x1) / width;
    var data = [];
    var i, j, x, y, row;

    //generate data to plot
    for(i = 0; i < width; i++) {
        x = x1 + i * xs;
        y = 1 + a2 / 2 * (3 * x * x - 1) + a4 / 8 * (35 * x * x * x * x - 30 * x * x + 3);
        row = [x];
        if(y.length > 0) {
            for(j = 0; j < y.length; j++) {
                row.push(y[j]);
            }
        } else {
            row.push(y);
        }
        data.push(row);
    }

    new Dygraph(graph, data,
        {
            xlabel: "cos &#952",
            ylabel: "W(&#952)",
            labels: ['Cos','W'],
            color: "red",
            strokeWidth: 3.0,
            valueRange: [0.0, 2.0]
        }
    );
};

function recalculate_l1values(){
    var j1 = document.getElementById("j1").value*1;
    var j2 = document.getElementById("j2").value*1;
    var j3 = document.getElementById("j3").value*1;
    var l1min, l1max, l1a, l1b, l1aradios, l1bradios, i;

    // l1 possibilities
    l1min = Math.abs(j2-j1);
    l1max = Math.abs(j2+j1);
    if (l1min==0){
        l1min = 1;
    }
    l1a = l1min;
    l1b = l1a+1;
    if (l1min==l1max){
        l1b = l1a;
    }

    l1aradios = "";
    for (i = l1min;i<=l1max;i++){
		    if (i==l1a)
		        l1aradios = l1aradios + "<input type=\"radio\" name=\"l1a\" value=" + i + " id=\"l1_" + i + "\" checked=\"checked\"><label for=\"l1_" + i + "\">" + i + "</label>";
		    else
		        l1aradios = l1aradios + "<input type=\"radio\" name=\"l1a\" value=" + i + " id=\"l1_" + i + "\"><label for=\"l1_" + i + "\">" + i + "</label>";
    }

    l1bradios = "";
    for (i = l1min;i<=l1max;i++){
		    if (i==l1b)
		        l1bradios = l1bradios + "<input type=\"radio\" name=\"l1b\" value=" + i + " id=\"l1_" + i + "\" checked=\"checked\"><label for=\"l1_" + i + "\">" + i + "</label>";
		    else
		        l1bradios = l1bradios + "<input type=\"radio\" name=\"l1b\" value=" + i + " id=\"l1_" + i + "\"><label for=\"l1_" + i + "\">" + i + "</label>";
    }
    document.getElementById("l1a_value").innerHTML = l1aradios;
    document.getElementById("l1b_value").innerHTML = l1bradios;
};

function recalculate_l2values(){
    var j2 = document.getElementById("j2").value*1;
    var j3 = document.getElementById("j3").value*1;
    var l2min, l2max, l2a, l2b, l2aradios, l2bradios, i;

    // l2 possibilities
    l2min = Math.abs(j3-j2);
    l2max = Math.abs(j3+j2);
    if (l2min==0){
        l2min = 1;
    }
    l2a = l2min;
    l2b = l2a+1;
    if (l2min==l2max){
        l2b = l2a;
    }

    l2aradios = "";
    for (i = l2min; i<=l2max; i++){
		    if (i==l2a)
		        l2aradios = l2aradios + "<input type=\"radio\" name=\"l2a\" value=" + i + " id=\"l2_" + i + "\" checked=\"checked\"><label for=\"l2_" + i + "\">" + i + "</label>";
		    else
		        l2aradios = l2aradios + "<input type=\"radio\" name=\"l2a\" value=" + i + " id=\"l2_" + i + "\"><label for=\"l2_" + i + "\">" + i + "</label>";
    }

    l2bradios = "";
      for (i = l2min; i<=l2max; i++){
          if (i==l2b)
		          l2bradios = l2bradios + "<input type=\"radio\" name=\"l2b\" value=" + i + " id=\"l2_" + i + "\" checked=\"checked\"><label for=\"l2_" + i + "\">" + i + "</label>";
		      else
		          l2bradios = l2bradios + "<input type=\"radio\" name=\"l2b\" value=" + i + " id=\"l2_" + i + "\"><label for=\"l2_" + i + "\">" + i + "</label>";
      }
      document.getElementById("l2a_value").innerHTML = l2aradios;
      document.getElementById("l2b_value").innerHTML = l2bradios;
};

function check_jvalues(){
        
    var j1 = document.getElementById("j1").value*1;
    var j2 = document.getElementById("j2").value*1;
    var j3 = document.getElementById("j3").value*1;

    document.getElementById("error1").innerHTML = "";

    if (j2==0 && (j1==0 || j3==0)){
        j2 = 1;
        document.getElementById("j2").value = 1*1;
        alert("No 0 to 0 transitions allowed. Setting J2 to 1");
    }
};

function recalculate(){
    var j1 = document.getElementById("j1").value*1;
    var j2 = document.getElementById("j2").value*1;
    var j3 = document.getElementById("j3").value*1;

    var l1a = $('input[name="l1a"]:checked').val()*1;
    var l1b = $('input[name="l1b"]:checked').val()*1;
    var l2a = $('input[name="l2a"]:checked').val()*1;
    var l2b = $('input[name="l2b"]:checked').val()*1;

    var d1 = $('#delta1-slider').attr('data-slider');
    var d2 = $('#delta2-slider').attr('data-slider');

    if (l1a==l1b){
        if (d1!=0){
            d1 = 0;
            $('.range-slider').foundation('delta1-slider', 'set_value', d1);
            alert("can't have mixing; only multipolarity selected is "+l1a);
        }
        $('#delta1-slider').addClass('disabled');
    } else {
        $('#delta1-slider').removeClass('disabled');
    }
    if (l2a==l2b){
		  if (d2!=0){
		    d2 = 0;
		    $('#delta2-slider').foundation('slider', 'set_value', d2);
		    alert("Can't have mixing; only multipolarity selected is "+l2a);
		  }
		  $('#delta2-slider').addClass('disabled');
    } else {
		  $('#delta2-slider').removeClass('disabled');
    }

    document.getElementById("a2").value = calculate_a2(j1,j2,j3,l1a,l1b,l2a,l2b,d1,d2);
    document.getElementById("a4").value = calculate_a4(j1,j2,j3,l1a,l1b,l2a,l2b,d1,d2);

    plot();
};

function calculate_a2(j1, j2, j3, l1a, l1b, l2a, l2b, delta1, delta2){
    return B(2,j2,j1,l1a,l1b,delta1)*A(2,j3,j2,l2a,l2b,delta2);
};

function calculate_a4(j1, j2, j3, l1a, l1b, l2a, l2b, delta1, delta2){
    return B(4,j2,j1,l1a,l1b,delta1)*A(4,j3,j2,l2a,l2b,delta2);
};

//------------------------------begin angular correlation functions--------------------//
function Factorial(value){
    var fac;
    if(value > 1){
        fac = value*Factorial(value-1);
    } else {
        fac = 1;
    }
    return fac;
}

function ClebschGordan(j1, m1, j2, m2, j, m){
    var term, cg, term1, sum, k;

    // Conditions check
    if( 2*j1 != Math.floor(2*j1) || 
        2*j2 !=   Math.floor(2*j2) || 
        2*j !=   Math.floor(2*j) || 
        2*m1 !=   Math.floor(2*m1) || 
        2*m2 !=   Math.floor(2*m2) || 
        2*m !=   Math.floor(2*m) ){

        //G4cout << "All arguments must be integers or half-integers." << G4endl;
        return 0;
    }

    if(m1 + m2 != m){
        //G4cout << "m1 + m2 must equal m." << G4endl;
        return 0;
    }

    if( j1 - m1 != Math.floor ( j1 - m1 ) ){
        //G4cout << "2*j1 and 2*m1 must have the same parity" << G4endl;
        return 0;
    }

    if( j2 - m2 != Math.floor ( j2 - m2 ) ){
        //G4cout << "2*j2 and 2*m2 must have the same parity" << G4endl;
        return 0;
    }

    if( j - m != Math.floor ( j - m ) ){
        //G4cout << "2*j and 2*m must have the same parity" << G4endl;
        return 0;
    }

    if(j > j1 + j2 || j < Math.abs(j1 - j2)){
        //G4cout << "j is out of bounds." << G4endl;
        return 0;
    }

    if(Math.abs(m1) > j1){
        //G4cout << "m1 is out of bounds." << G4endl;
        return 0;
    }

    if(Math.abs(m2) > j2){
        //G4cout << "m2 is out of bounds." << G4endl;
        return 0;
    }

    if(Math.abs(m) > j){
        //warning('m is out of bounds." << G4endl;
        return 0 ;
    }

    term1 = Math.pow((((2*j+1)/Factorial(j1+j2+j+1))*Factorial(j2+j-j1)*Factorial(j+j1-j2)*Factorial(j1+j2-j)*Factorial(j1+m1)*Factorial(j1-m1)*Factorial(j2+m2)*Factorial(j2-m2)*Factorial(j+m)*Factorial(j-m)),(0.5));
    sum = 0;
    
    for(k = 0 ; k <= 99 ; k++ ){
        if( (j1+j2-j-k < 0) || (j-j1-m2+k < 0) || (j-j2+m1+k < 0) || (j1-m1-k < 0) || (j2+m2-k < 0) ){}
        else{
            term = Factorial(j1+j2-j-k)*Factorial(j-j1-m2+k)*Factorial(j-j2+m1+k)*Factorial(j1-m1-k)*Factorial(j2+m2-k)*Factorial(k);
            if((k%2) == 1){
                term = -1*term;
            }
            sum = sum + 1.0/term;
        }
    }
    
    cg = term1*sum;
    return cg;
    // Reference: An Effective Algorithm for Calculation of the C.G.
    // Coefficients Liang Zuo, et. al.
    // J. Appl. Cryst. (1993). 26, 302-304
};

function Wigner3j(j1, j2, j3, m1, m2, m3){
    var out;

    // Conditions check
    if( 2*j1 != Math.floor(2*j1) || 
        2*j2 != Math.floor(2*j2) || 
        2*j3 != Math.floor(2*j3) || 
        2*m1 != Math.floor(2*m1) || 
        2*m2 != Math.floor(2*m2) || 
        2*m3 != Math.floor(2*m3) ){
        // G4cout << "All arguments must be integers or half-integers." << G4endl;
        return 0;
    }

    if(m1 + m2 + m3 != 0){
        //G4cout << "m1 + m2 + m3 must equal zero." << G4endl;
        return 0;
    }

    if( j1 + j2 + j3 !=   Math.floor(j1 + j2 + j3) ){
        //G4cout << "2*j1 and 2*m1 must have the same parity" << G4endl;
        return 0;
    }

    if(j3 > j1 + j2 || j3 < Math.abs(j1 - j2)){
        //G4cout << "j3 is out of bounds." << G4endl;
        return 0;
    }

    if(Math.abs(m1) > j1){
        //G4cout << "m1 is out of bounds." << G4endl;
        return 0;
    }

    if(Math.abs(m2) > j2){
        //G4cout << "m2 is out of bounds." << G4endl;
        return 0;
    }

    if(Math.abs(m3) > j3){
        return 0;
    }

    out = (Math.pow((-1),(j1-j2-m3)))/(Math.pow((2*j3+1),(1.0/2.0)))*ClebschGordan(j1,m1,j2,m2,j3,-1*m3);
        return out;
};

function Wigner6j(J1, J2, J3, J4, J5, J6){
    var j1 = J1;
    var j2 = J2;
    var j12 = J3;
    var j3 = J4;
    var j = J5;
    var j23 = J6;
    var sum = 0;

    // Conditions check
    if(J3 > J1 + J2 || J3 < Math.abs(J1 - J2)){
        //G4cout << "first J3 triange condition not satisfied. J3 > J1 + J2 || J3 < Math.abs(J1 - J2)" << G4endl;
        return 0;
    }

    if(J3 > J4 + J5 || J3 < Math.abs(J4 - J5)){
        //G4cout << "second J3 triange condition not satisfied. J3 > J4 + J5 || J3 < Math.abs(J4 - J5)" << G4endl;
        return 0;
    }

    if(J6 > J2 + J4 || J6 < Math.abs(J2 - J4)){
        //G4cout << "first J6 triange condition not satisfied. J6 > J2 + J4 || J6 < Math.abs(J2 - J4)" << G4endl;
        return 0;
    }

    if(J6 > J1 + J5 || J6 < Math.abs(J1 - J5)){
        //G4cout << "second J6 triange condition not satisfied. J6 > J1 + J5 || J6 < Math.abs(J1 - J5)" << G4endl;
        return 0;
    }
    
    for(var m1 = -j1 ; m1 <= j1 ; m1++ ){
        for(var m2 = -j2 ; m2 <= j2 ; m2++ ){
            for(var m3 = -j3 ; m3 <= j3 ; m3++ ){
                for(var m12 = -j12 ; m12 <= j12 ; m12++ ){
                    for(var m23 = -j23 ; m23 <= j23 ; m23++ ){
                        for(var m = -j ; m <= j ; m++ ){
                            sum = sum + Math.pow((-1),(j3+j+j23-m3-m-m23))*Wigner3j(j1,j2,j12,m1,m2,m12)*Wigner3j(j1,j,j23,m1,-m,m23)*Wigner3j(j3,j2,j23,m3,m2,-m23)*Wigner3j(j3,j,j12,-m3,m,m12);
                        }
                    }
                }
            }
        }
    }
    return sum;
};

function RacahW(a, b, c, d, e, f){
    return Math.pow((-1),(a+b+d+c))*Wigner6j(a,b,e,d,c,f);
};

function F(k, jf, L1, L2, ji){
    var W;
    var CG = ClebschGordan(L1,1,L2,-1,k,0);

    if(CG == 0){
        return 0;
    }
    W = RacahW(ji,ji,L1,L2,k,jf);
    if(W == 0){
        return 0;
    }
    return Math.pow((-1),(jf-ji-1))*(Math.pow((2*L1+1)*(2*L2+1)*(2*ji+1),(1.0/2.0)))*CG*W;
    // Reference: Tables of coefficients for angular distribution of gamma rays from aligned nuclei
    // T. Yamazaki. Nuclear Data A, 3(1):1?23, 1967.
};

function A(k, ji, jf, L1, L2, delta){
    return (1/(1+Math.pow(delta,2)))*(F(k,ji,L1,L1,jf)+2*delta*F(k,ji,L1,L2,jf)+delta*delta*F(k,ji,L2,L2,jf));
};

function B(k, ji, jf, L1, L2, delta){
      return (1/(1+Math.pow(delta,2)))*(  F(k,jf,L1,L1,ji)+(Math.pow((-1),((L1+L2))))*2*delta*F(k,jf,L1,L2,ji)+delta*delta*F(k,jf,L2,L2,ji) );
};

function evenA(){
    
    var select, option, i, spin;

    for(spin = 1; spin<4; spin++){
        select = document.getElementById('j'+spin);
        select.innerHTML = '';
        for (i = 0;i<10;i++){
                option = document.createElement('option');
                option.setAttribute('value', i);
                option.innerHTML = i;
                select.appendChild(option);
        }
    }

    document.getElementById("j1").value = 4;
    document.getElementById("j2").value = 2;
    document.getElementById("j3").value = 0;
};

function oddA(){
    var select, option, i, spin;

    for(spin = 1; spin<4; spin++){
        select = document.getElementById('j'+spin);
        select.innerHTML = '';
        for (i = 0;i<10;i++){
                value = (2*i+1)/2;
                numerator = 2*i+1;
                
                option = document.createElement('option');
                option.setAttribute('value', value);
                option.innerHTML = numerator + '/2';
                select.appendChild(option);                
        }  
    }

    document.getElementById("j1").value = 2.5;
    document.getElementById("j2").value = 1.5;
    document.getElementById("j3").value = 0.5;
};

function plot(setup){
    var graph = document.getElementById("graph_div"),
        width = 1000,
        x1 = -1,
        x2 = 1,
        a2 = parseFloat(document.getElementById("a2").value),
        a4 = parseFloat(document.getElementById("a4").value),
        xs = 1.0 * (x2 - x1) / width,
        data = [],
        plotWidth = document.getElementById('plotCol').offsetWidth,
        plotHeight = plotWidth*3/4,
        i, j, x, y, row;

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

    //dygraphs fixes div size on paint, remove to allow resize
    document.getElementById('graph_div').setAttribute('style', '');

    //first time setup the dygraphs
    if(setup){
        dataStore.plot = new Dygraph(graph, [[0,0],[0,0]],
            {
                xlabel: "cos &#952",
                ylabel: "W(&#952)",
                labels: ['Cos','W'],
                color: "red",
                strokeWidth: 3.0,
                valueRange: [0.0, 2.0],
                width: plotWidth,
                height: plotHeight
            }
        );
    } else {
        dataStore.plot.updateOptions( { 'file': data });
    }
};

function plot2D(n){
    //regenerate the plot for a2 (n=2) or a4 (n=4)
    var data = [
            {
                x: dataStore.x,
                y: dataStore.y,
                z: (n==2) ? dataStore.a2 : dataStore.a4,
                type: 'contour',
                name: (n==2) ? 'a2': 'a4', 
                hoverinfo:"x+y+z",
                colorscale: 'Viridis'
            }
        ],
        dim = document.getElementById('a'+n+'Wrap').offsetWidth,
        layout = {
            title: (n==2) ? 'a2' : 'a4',
            xaxis:{
                title: 'L<sub>1a</sub> / L<sub>1b</sub> mixing'
            },
            yaxis:{
                title: 'L<sub>2a</sub> / L<sub>2b</sub> mixing'
            },
            autosize: false,
            width: dim,
            height: dim
        }


    Plotly.newPlot('a'+n+'Plot', data, layout);
}

function plot_a(n, missingL1, missingL2){
    //regenerate the 1D plot for at (n==2) or a4 (n==4)
    //leave the appropriate mixing ratio unconstrained depending on who's not missing (should never be both missing)

    var data,
        dim = document.getElementById('a'+n+'Wrap').offsetWidth,
        layout = {
            autosize: false,
            width: dim,
            height: dim,
            xaxis:{
                title: missingL1 ? 'L<sub>2a</sub> / L<sub>2b</sub> mixing' : 'L<sub>1a</sub> / L<sub>1b</sub> mixing'
            },
            yaxis:{
                title: (n==2) ? 'a2' : 'a4'
            },
        }, i,
        zero = 1/(2/dataStore.steps)

    //slice this 2D data along 0 for the missing mixing ratio 
    if(missingL1){
        data = [];
        for(i=0; i<dataStore.steps; i++){
            data.push( (n==2) ? dataStore.a2[i][zero] : dataStore.a4[i][zero] );
        }
    } else if(missingL2){
        data = (n==2) ? dataStore.a2[zero] : dataStore.a4[zero]
    }

    //construct the plotly data object
    data = [
        {
            x: dataStore.x,
            y: data,
            name: (n==2) ? 'a2': 'a4',
            mode: 'markers',
            type: 'scatter'
        }
    ]
    
    Plotly.newPlot('a'+n+'Plot', data, layout);    
}

function recalculate_L(transition){
    // transition == 1 -> first transition; 2 -> second transition.

    var jOrig = parseFloat(document.getElementById("j" + transition).value),
        jFinal = parseFloat(document.getElementById("j" + (transition+1)).value),
        lMin, lMax, lA, lB, i, j,
        radio, label, chosenMomenta,
        momenta = ['l'+transition+'a', 'l'+transition+'b'],
        momenta_options = [document.getElementById('l'+transition+'a_value'), document.getElementById('l'+transition+'b_value')];

    lMin = Math.abs(jFinal-jOrig);
    lMax = Math.abs(jFinal+jOrig);
    if (lMin==0){
        lMin = 1;
    }
    lA = lMin;
    lB = lA+1;
    if (lMin==lMax){
        lB = lA;
    }
    chosenMomenta = [lA, lB];

    for(j=0; j<momenta.length; j++){
        momenta_options[j].innerHTML = '';

        label = document.createElement('label');
        label.innerHTML = "L<sub>" + momenta[j].slice(1) + "</sub>";
        momenta_options[j].appendChild(label);

        for (i = lMin; i<=lMax; i++){
            radio = document.createElement('input');
            radio.setAttribute('type', 'radio');
            radio.setAttribute('name', momenta[j]);
            radio.setAttribute('value', i);
            radio.setAttribute('id', momenta[j] + '_' + i);
            radio.onchange = recalculate;
            if(i == chosenMomenta[j])
                radio.setAttribute('checked', true);
            momenta_options[j].appendChild(radio);

            label = document.createElement('label');
            label.setAttribute('for', momenta[j] + '_' + i);
            label.innerHTML = i;
            momenta_options[j].appendChild(label);
        }
    } 
};

function check_jvalues(){
        
    var j1 = parseFloat(document.getElementById("j1").value);
        j2 = parseFloat(document.getElementById("j2").value),
        j3 = parseFloat(document.getElementById("j3").value),

    document.getElementById("error1").innerHTML = "";

    if (j2==0 && (j1==0 || j3==0)){
        j2 = 1;
        document.getElementById("j2").value = 1;
        alert("No 0 to 0 transitions allowed. Setting J2 to 1");
    }
};

function recalculate(){
    
    var j1 = parseFloat(document.getElementById("j1").value),
        j2 = parseFloat(document.getElementById("j2").value),
        j3 = parseFloat(document.getElementById("j3").value),

        l1a = parseFloat($('input[name="l1a"]:checked').val()),
        l1b = parseFloat($('input[name="l1b"]:checked').val()),
        l2a = parseFloat($('input[name="l2a"]:checked').val()),
        l2b = parseFloat($('input[name="l2b"]:checked').val()),

        d1 = parseFloat($('#mix1').val()),
        d2 = parseFloat($('#mix2').val()),
        i, j, row,
        noL1mix = false, 
        noL2mix = false;

    if (l1a==l1b){
        noL1mix = true;
        if (d1!=0){
            d1 = 0;
            $('#mix1').val(d2);
            $('#delta1-slider').val(d2);
            alert("can't have mixing; only multipolarity selected is "+l1a);
        }
        $('#delta1-slider').attr('disabled', 'true');
    } else {
        $('#delta1-slider').removeAttr('disabled');
    }
    if (l2a==l2b){
        noL2mix = true;
		if (d2!=0){
		    d2 = 0;
		    $('#mix2').val(d2);
            $('#delta2-slider').val(d2);
		    alert("Can't have mixing; only multipolarity selected is "+l2a);
		}
		$('#delta2-slider').attr('disabled', 'true');
    } else {
		$('#delta2-slider').removeAttr('disabled');
    }

    document.getElementById("a2").value = calculate_a2(j1,j2,j3,l1a,l1b,l2a,l2b,d1,d2);
    document.getElementById("a4").value = calculate_a4(j1,j2,j3,l1a,l1b,l2a,l2b,d1,d2);

    plot();

    document.getElementById('customAwarning').classList.add('hidden');

    //a2 and a4 plots
    //generate data
    dataStore.a2 = [];
    dataStore.a4 = [];
    for(i=0; i<dataStore.steps; i++){
        dataStore.a2[i] = []
        dataStore.a4[i] = []
        for(j=0; j<dataStore.steps; j++){
            dataStore.a2[i][j] = dataStore.A2[i]*dataStore.B2[j];
            dataStore.a4[i][j] = dataStore.A4[i]*dataStore.B4[j];
        }
    }

    if (noL1mix && noL2mix){
        document.getElementById('a2Plot').innerHTML = '';
        document.getElementById('a4Plot').innerHTML = '';
    }
    else if(noL1mix || noL2mix){
        plot_a(2, noL1mix, noL2mix);
        plot_a(4, noL1mix, noL2mix);
    } else {
        plot2D(2);
        plot2D(4);
    }
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

    if(dataStore.cache.factorial[value]){
        return dataStore.cache.factorial[value];
    } else {
        if(value > 1){
            fac = value*Factorial(value-1);
        } else {
            fac = 1;
        }
        dataStore.cache.factorial[value] = fac;

        return fac;
    }
}

function ClebschGordan(j1, m1, j2, m2, j, m){
    var term, cg, term1, sum, k

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
        if( (j1+j2-j-k < 0) || (j1-m1-k < 0) || (j2+m2-k < 0) )
            //no further terms will contribute to sum, exit loop
            break
        else if( (j-j1-m2+k < 0) || (j-j2+m1+k < 0)  )
            //jump ahead to next term that will contribute
            k = Math.max(Math.max(j-j1-m2, j-j2+m1) - 1, k);
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
    // if( 2*j1 != Math.floor(2*j1) || 
    //     2*j2 != Math.floor(2*j2) || 
    //     2*j3 != Math.floor(2*j3) || 
    //     2*m1 != Math.floor(2*m1) || 
    //     2*m2 != Math.floor(2*m2) || 
    //     2*m3 != Math.floor(2*m3) ){
    //     // G4cout << "All arguments must be integers or half-integers." << G4endl;
    //     return 0;
    // }

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
        j2 = J2,
        j12 = J3,
        j3 = J4,
        j = J5,
        j23 = J6,
        sum = 0;

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
    var f1 = F(k,ji,L1,L1,jf),
        f2 = F(k,ji,L1,L2,jf),
        f3 = F(k,ji,L2,L2,jf);

    tabulateA(k, f1,f2,f3);

    return (1/(1+Math.pow(delta,2)))*(f1+2*delta*f2+delta*delta*f3);
};

function tabulateA(k, f1, f2, f3){
    //given precomputed values of F, reconstruct the table of A values for the currently selected momenta, across a range of mixing ratios.
    var i, delta;

    if(k==2)
        dataStore.A2 = [];
    else if(k==4)
        dataStore.A4 = []
    for(i=0; i<dataStore.steps; i++){
        delta = -1 + 2*i/dataStore.steps;
        if(k==2)
            dataStore.A2.push( (1/(1+Math.pow(delta,2)))*(f1+2*delta*f2+delta*delta*f3) );
        else if(k==4)
            dataStore.A4.push( (1/(1+Math.pow(delta,2)))*(f1+2*delta*f2+delta*delta*f3) );
    }
}

function B(k, ji, jf, L1, L2, delta){
    var f1 = F(k,jf,L1,L1,ji),
        f2 = F(k,jf,L1,L2,ji),
        f3 = F(k,jf,L2,L2,ji)

    tabulateB(k, f1,f2,f3,L1,L2);

    return (1/(1+Math.pow(delta,2)))*(f1+(Math.pow((-1),((L1+L2))))*2*delta*f2+delta*delta*f3);
};

function tabulateB(k, f1, f2, f3, L1, L2){
    //given precomputed values of F, reconstruct the table of B values for the currently selected momenta, across a range of mixing ratios.
    var i, delta;

    if(k==2)
        dataStore.B2 = [];
    else if(k==4)
        dataStore.B4 = [];
    for(i=0; i<dataStore.steps; i++){
        delta = -1 + 2*i/dataStore.steps;
        if(k==2)
            dataStore.B2.push( (1/(1+Math.pow(delta,2)))*(f1+(Math.pow((-1),((L1+L2))))*2*delta*f2+delta*delta*f3) );
        else if (k==4)
            dataStore.B4.push( (1/(1+Math.pow(delta,2)))*(f1+(Math.pow((-1),((L1+L2))))*2*delta*f2+delta*delta*f3) );
    }
}

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

/////////////////
// ui helpers
/////////////////

function syncElements(source, dest){
    //source: string; id of element to read value from
    //dest: string; id of element to write value to
    //onchange callback to set the value of another element to that of this one.

    var val = document.getElementById(source).value;
    document.getElementById(dest).value = val;
}

Release History:

version | DOI
--------|------
1.0.0   | 

# Angular Correlation Calculator

The Angular Correlation Calculator takes spins, multipolarities, and mixing ratios for a gamma-gamma cascade, and calculates the angular correlation coefficients and plots the correlation. Also, plots of a2 and a4 as functions of one or both mixing ratios are generated as appropriate, as well as a parametric plot of a2 vs. a4 as a function of mixing ratio, when only one mixing ratio remains unconstrained.

## Dependencies & Setup

The Angular Correlation Calculator runs 100% client side; simply open `index.html` in the latest Firefox or Chrome locally, or serve from any static page server.

This project uses [Dygraphs](http://dygraphs.com/) and [Plotly](https://plot.ly/) for plotting, [Twitter Bootstrap](http://getbootstrap.com/) for layout, and [QUnit](http://qunitjs.com/) for unit testing.
 
## Programmatic Logic

A few key points in how this calculator is put together:

 - All global variables are namespaces on a `dataStore` object, created at page load; **this is the appropriate place to add other global variables.**
 - All calculations and state changes are triggered in response to user interactions with the UI panel. Have a look at the `oninput` and/or `onchange` properties inline on the HTML of these controls to see what functions get called in response to changing them. In general, angular correlations and available spin / multipolarity options are recalculated, and plots are regenerated bases on the new state of the UI. 

## Contributing

Contributions are very welcome! If you have an idea, question or comment, please open an issue. If you would like to make a change to this project, please follow these steps:
 - start by opening an issue or empty PR to discuss your ideas
 - please limit individual PRs to less than 500 lines (Why? See figure 1 [here](https://smartbear.com/SmartBear/media/pdfs/11_Best_Practices_for_Peer_Code_Review.pdf)).
 - please encapsulate all new behavior wherever possible in functions of 50 lines or less each.
 - please include unit tests for all new functions wherever possible, to demonstrate correct behavior.

## Citation & Deployment

If you use a result from this project, **be sure to site it using the correct DOI**. This will allow you to go back and reproduce your results later, with the same software version you used originally. To find the correct DOI, look in the footer of the app.

If you push changes to this project onto GRIFFIN's live toolkit, **be sure to update the DOI in the footer and in the table at the top of this file**. To get a new DOI, simply [make a new release via GitHub](https://help.github.com/articles/creating-releases/), then [visit Zenodo](https://zenodo.org/account/settings/github/), sign in with your GitHub credentials, and find this project in the list on that page; clicking on the badge will give you a bunch of options to cut and paste into the appropriate places. Add the markdown one to this document, and the HTML one to the footer.






# AngularCorrelationUtility
Web page of useful tools for gamma-gamma angular correlations
* Given spins, multipolarities, and mixing ratios for a gamma-gamma cascade, this page calculates the angular correlation coefficients and displays the correlation.
* Uses functions from Geant4GammaGammaAngularCorrelations10.01.
* Also uses the dygraphs Javascript library and the Foundation Javascript library and CSS formatting.

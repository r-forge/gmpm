
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<h2>Generalized Multilevel Permutation Models</h2>

<p>The gmpm package provides a comprehensive framework for performing permutation tests using regression on multilevel experimental data.  GMPMs can be used to analyze categorical, count, and continuous data on its natural scale.</p>

<!-- end of project description -->

<p>The package <i>gmpm</i> can be installed from the R command line using the following syntax:</p>

<p><strong><code>install.packages("gmpm", repos="http://R-Forge.R-project.org")</code></strong></p>

  <p>To get started with <i>gmpm</i>, please see the package help documentation for the function <i>gmp</i>.  If you are interested in using GMPM to analyze visual-world eyetracking data, please see the example <i>kb07</i> in the package (from the command line in R, type <code>?kb07</code>).</p>

  <p>New features will be added soon.  If you want to receive email updates, please subscribe to the <a href="https://lists.r-forge.r-project.org/cgi-bin/mailman/listinfo/gmpm-news">'news' mailing list</a>.</p>

  <p>Please see the <a href="https://r-forge.r-project.org/forum/?group_id=512">forums</a> for news/discussion.  If you encounter problems, please submit a report in the "bugs" forum.</p>

  <p> You can find the <strong>project summary page</strong> <a href="http://r-forge.r-project.org/projects/gmpm"><strong>here</strong></a>. </p>

</body>
</html>

<?php

/**
 * @file
 *
 */

/**
 * generate form for go enrichment analysis.
 */
function pwy_tool_enrich_form($form, $form_state) {

  //load default value from RNASeq 
  $default_lists = '';
  if(module_exists('tripal_rnaseq')) {
    if (isset($_SESSION['tripal_rnaseq_analysis']['genes_str'])) {
      $default_lists = $_SESSION['tripal_rnaseq_analysis']['genes_str'];
    }
    if (isset($_SESSION['tripal_rnaseq_analysis']['organism_id'])) {
      $organism_id = $_SESSION['tripal_rnaseq_analysis']['organism_id'];
      $organism = chado_generate_var('organism', array('organism_id'=>$organism_id));
    }
  }

  //* choose pwy dataset
  //*   An array where the nid is the key and the value is array
  //*   array: 
  //*      1 the human-readable name of the option
  //*      2 path
  //*      3 link for gene id
  //*      4 link for pwy id
  //*      5 example id
  $pwy_db = get_pwy_db_options();

  $id_link = array();
  $db_link = array();
  $db_expid = array();
  $options = array(
    'SELECT_DB' => 'Select a Dataset',
  );

  $default_opt = 'SELECT_DB';

  foreach ($pwy_db as $nid => $pwy_array) {
    $options[$pwy_array[1]]  = $pwy_array[0];
    if ($pwy_array[0] == $organism->common_name) {
      $default_opt = $pwy_array[1];
    }
    $id_link[$pwy_array[1]]  = $pwy_array[2];
    $db_link[$pwy_array[1]]  = $pwy_array[3];
    $db_expid[$pwy_array[1]] = $pwy_array[4];
  }

  // pass these array for ajax, update type, link, and example id
  $form['#id_link']  = $id_link;
  $form['#db_link']  = $db_link;
  $form['#db_expid'] = $db_expid;

  $form['SELECT_PWY_DATASET'] = array(
    '#type' => 'select',
    '#title' => t('Dataset:'),
    '#options' => $options,
    '#default_value' => $default_opt,
    '#ajax' => array(
      'event' => 'change',
      'callback' => 'ajax_pwy_tool_enrich_example_id_callback',
      'method'   => 'replace',
      'wrapper'  => 'id-textarea',
      'effect'   => 'fade',
    ),
  );

  // NUCLEOTIDE QUERY
  //.........................
  $form['query'] = array(
    '#type' => 'fieldset',
    '#title' => t('Enter the list of IDs.'),
    '#collapsible' => TRUE,
    '#collapsed' => FALSE,
  );

  // Checkbox to show an example.
  $form['query']['example_id'] = array(
    '#type' => 'checkbox',
    '#title' => t('Show an Example ID'),
    '#prefix' => '<span style="float: right;">',
    '#suffix' => '</span>',
    '#ajax' => array(
      'callback' => 'ajax_pwy_tool_enrich_example_id_callback',
      'wrapper'  => 'id-textarea',
      'method'   => 'replace',
      'effect'   => 'fade',
    ),
  );

  // Textfield for submitting a ID list for query
  $form['query']['ID'] = array(
    '#type' => 'textarea',
    '#title' => t('Enter gene IDs'),
    '#description'=>t('Enter gene IDs in the text area.'),
    '#default_value' => $default_lists,
    '#prefix' => '<div id="id-textarea">',
    '#suffix' => '</div>',
  );

  // Upload a file as an alternative to enter a query ID
  $form['#attributes']['enctype'] = 'multipart/form-data';
  $form['query']['UPLOAD'] = array(
    '#title' => 'Or upload the file containing the list of gene IDs:  ',
    '#type' => 'managed_file',
    '#description' => t('The file should be a plain-text
  (.txt) file. Max size: %max_size. <strong>Don\'t forget to press the Upload button before
  attempting to submit your ID.</strong>',
      array(
        '%max_size' => round(file_upload_max_size() / 1024 / 1024,1) . 'MB'
      )
    ),
    '#upload_validators' => array(
      'file_validate_extensions' => array('txt'),
      'file_validate_size' => array(file_upload_max_size()),
    ),
  );

  $form['cutoff_pvalue'] = array(
    '#type' => 'textfield',
    '#title' => t('Enter a cutoff P-Value:'),
    '#default_value' => '0.05',
    '#maxlength' => 20,
  );

  // Submit
  //.........................
  $form['submit'] = array(
    '#type' => 'submit',
    '#default_value' => ' Submit ',
  );

  $form['#prefix'] = '<div class="row"> <div class="col-md-8 col-md-offset-2">';
  $form['#suffix'] = '</div> </div>';

  return $form;
}

/**
 * Validate the user options submitted via the above form.
 *
 * @see pwy_tool_enrhich_form().
 */
function pwy_tool_enrich_form_validate($form, &$form_state) {

  // Validate Query
  //----------------
  // @todo: We are currently not validating uploaded files are valid FASTA.
  // First check to see if we have an upload & if so then validate it.
  $file = file_load($form_state['values']['UPLOAD']);

  // If the $file is populated then this is a newly uploaded, temporary file.
  if (is_object($file)) {
    $form_state['qFlag'] = 'upQuery';
    $form_state['upQuery_path'] = drupal_realpath($file->uri);
  }
  // Otherwise there was no file uploaded.
  // Check if there was a query sequence entered in the texfield.
  elseif (!empty($form_state['input']['ID'])) {
    // Check to ensure that the query sequence entered is valid FASTA.
    // if (validate_fasta_sequence($type, $form_state['input']['ID'])){
    //  form_set_error('query', t('You need to provide a valid FASTA sequence '
    //  . 'for the query. It should contain a FASTA header/definition line '
    //  . 'followed by %molecule-type. For more information see the '
    //  . '<a href="@url" target="_blank">NCBI FASTA Specification</a>.',
    //    array(
    //      '%molecule-type' => $molecule_type,
    //      '@url' => 'http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml'
    //  )));
    //}
    //else {
      $form_state['qFlag'] = 'seqQuery';
    //}
  }
  // Otherwise they didn't enter a query!!
  else {
    form_set_error('query', t('No gene IDs given. Enter IDs in the box provided '
    . 'or upload a plain text file.'));
  }

  // Validate Dataset
  if (!empty($form_state['values']['SELECT_PWY_DATASET'])) {
  }
  // Otherwise they didn't select a database!!
  else {
    form_set_error('DB', t('No dataset selected. Choose a dataset from the list.'));
  }

  // Validate ID list
  if (!empty($form_state['values']['ID']) && $form_state['values']['ID'] == 'unknown dataset type') {
    form_set_error('query', t('No gene IDs given.')); 
  } 

  //dpm($form_state['values']['SELECT_PWY_DATASET']);
  //form_set_error('query', t('No gene IDs given.'));

}//go_tool_enrich_form_validate


/**
 * Process the user options submitted.
 *
 * @see pwy_tool_enrich_form().
 */
function pwy_tool_enrich_form_submit($form, &$form_state) {

  $error = FALSE;

  // First set defaults.
  $pwyenrichjob = array(
    'job_id' => NULL,
    'pwy_file' => (isset($form_state['values']['SELECT_PWY_DATASET'])) ? $form_state['values']['SELECT_PWY_DATASET'] : NULL,
    'gene_list' => NULL,
    'result_file' => NULL,
    'options' => serialize(array())
  );

  // QUERY
  //-------------------------
  // BLAST can only take the query as a file;
  // therefore, if it was submitted via the textfield we need to create a file containing
  // the submitted sequence.
  if (isset($form_state['qFlag'])) {
    if ($form_state['qFlag'] == 'seqQuery') {
      $seq_content = $form_state['values']['ID'];
      $pwyenrichjob['gene_list'] = '/tmp/' . date('YMd_His') . '_pwy_gene_list.txt';
      file_put_contents ($pwyenrichjob['gene_list'], $seq_content);
    }
    elseif ($form_state['qFlag'] == 'upQuery') {
      $pwyenrichjob['gene_list'] = $form_state['upQuery_path'];
    }
  }

  // pathway annotation dataset
  $selected_db = $form_state['values']['SELECT_PWY_DATASET'];
  $pwyenrichjob['pwy_file'] = $selected_db; //GAF_dataset_path_full;


  

  // Determine the path to the pwy with extension.
  if (!is_readable($selected_db)) {
    $error = TRUE;
    // Reports error in admin>>reports>>Recent Log Messages
    tripal_report_error(
      'pwy_tool',
      TRIPAL_ERROR,
      "The PWY dataset is unaccessible."
    );
    // report error in current page
    $msg = "The pwy dataset '$selected_db' is unaccessible. ";
    $msg .= "Please contact the site administrator.";
    drupal_set_message($msg, 'error');
  }

  // OPTIONS
  $options = array();
  if (isset($form['#id_link'][$selected_db]) && isset($form['#db_link'][$selected_db])) {
    $options = array(
      't' => $form['#id_link'][$selected_db],
      'x' => $form['#db_link'][$selected_db],
      'p' => $form_state['values']['cutoff_pvalue'],
    );

  } else {
    tripal_report_error(
      'pwy_tool',
      TRIPAL_ERROR,
      "The PWY type and link are undetermined."
    );
  }

  $pwyenrichjob['options'] = serialize($options);
  $pwyenrich_program = $goenrich_program = DRUPAL_ROOT . DIRECTORY_SEPARATOR . drupal_get_path('module', 'pwy_tool') . DIRECTORY_SEPARATOR . 'pwy_tool.pl';

  // SUBMIT JOB TO TRIPAL
  //-------------------------
  if (!$error) {
    // pwy dataset exists.
    global $user;

    // set output dir and file
    $output_dir = variable_get('file_public_path', conf_path() . '/files')
      . DIRECTORY_SEPARATOR . 'tripal' . DIRECTORY_SEPARATOR . 'tripal_pathway';
    $output_file = $output_dir . DIRECTORY_SEPARATOR . date('YMd_His') . '.pwyenrich.txt';

    $job_args = array(
      'program' => $pwyenrich_program,
      'gene_list' => $pwyenrichjob['gene_list'],
      'pwy_file' => $pwyenrichjob['pwy_file'],
      'output_file' => $output_file,
      'options' => $options
    );

    $job_id = tripal_add_job(
      t('Pathway enrichment for @pwy_file', array('@pwy_file' => $pwyenrichjob['pwy_file'])),
      'tripal_pathway',
      'run_pwyenrich_tripal_job',
      $job_args,
      $user->uid
    );

    $pwyenrichjob['result_file'] = $output_file;
    $pwyenrichjob['job_id'] = $job_id;

    //Encode the job_id
    // the encode id is ordered in tripal jobs
    $job_encode_id = pwy_tool_make_secret($job_id);

    // RECENT JOBS
    //-------------------------
    if (!isset($_SESSION['pwyenrich_jobs'])) {
      $_SESSION['pwyenrich_jobs'] = array();
    }
    $_SESSION['pwyenrich_jobs'][] = $job_encode_id;

    // NOTE: Originally there was a call to tripal_launch_jobs() here. That should
    // NEVER be done since it runs possibly long jobs in the page load causing time-out
    // issues. If you do not want to run tripal jobs manually, look into installing
    // Tripal daemon which will run jobs as they're submitted or set up a cron job to
    // launch the tripal jobs on a specified schedule.

    // Redirect to the BLAST results page
    drupal_goto("pwyenrich/report/$job_encode_id");
  }
}

/**
 * AJAX: Replace the ID textarea with one containing an example.
 */
function ajax_pwy_tool_enrich_example_id_callback($form, $form_state) {
  $dataset = $form_state['values']['SELECT_PWY_DATASET'];
  $example_ids = 'unknown pwy type';

  if (isset($form['#db_expid'][$dataset])) {
    $example_ids = $form['#db_expid'][$dataset];
  }

  // If the Show Example checkbox is true then put the example in the textfield
  if (!isset($_SESSION['tripal_rnaseq_analysis']['genes_str'])) {
    if ($form_state['values']['example_id']) {
      $form['query']['ID']['#value'] = $example_ids;
    }
    else {
      $form['query']['ID']['#value'] = '';
    }
  }

  return $form['query']['ID'];
}


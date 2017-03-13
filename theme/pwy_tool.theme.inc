<?php

/**
 * @file
 * This file contains functions specifically related to theme-ing the BLAST module
 */


/**
 * Preprocess function for the show_pwyenrich_report.tpl.php
 *
 * Use this function to prepare variables for use in the template,
 * as well as to add css/javascript needed.
 *
 * @param $vars
 *   The variables currently available to the template.
 *
 * after preprocess: 
 * the var tripal_job, html_path, html_info could be directly use in them: show_goenrich_report
 *   theme show_goenrich_report point to template: goenrich_report
 *   so, these var is directly used in goenrich_report template
 */
function pwy_tool_preprocess_show_pwyenrich_report(&$vars) {
  global $base_url;

  // Add CSS and Javascript files
  $path = drupal_get_path('module', 'pwy_tool');
  // drupal_add_css($path . '/theme/css/your_css.css');
  // drupal_add_js('https://code.jquery.com/jquery-1.12.4.min.js');

  $tripal_job = tripal_get_job($vars['job_id']);
  $vars['tripal_job'] = $tripal_job;
  $tripal_args = unserialize($tripal_job->arguments);
  $vars['html_path'] = DRUPAL_ROOT . DIRECTORY_SEPARATOR . $tripal_args['output_file'];
  $vars['html_info'] = 'No result found';

  # load the html info to var
  if (is_readable($vars['html_path'])) {
   
    $result_table = array();
    $fp = fopen($vars['html_path'],'r');
    if (($headers = fgetcsv($fp, 0, "\t")) !== FALSE) {
      if ($headers) {
        while (($line = fgetcsv($fp, 0, "\t")) !== FALSE) {
          if ($line) {
            if (sizeof($line)==sizeof($headers)) {
              $result_table[] = array_combine($headers,$line);
            } 
            else {
              if ($line[0] == '#idlink') {  
                $vars['idlink'] = $line[1];
                if (preg_match('/[site_url]/', $line[1])) {
                  $vars['idlink'] = substr($base_url, 7) . substr($line[1], 10);
                }
              }
              if ($line[0] == '#dblink') {
                $vars['dblink'] = $line[1];
              }
            }
          }
        }
      }
    }
    fclose($fp);

    $vars['result_table'] = $result_table;
  }
}

/**
 * Implements hook_theme_registery_alter().
 *
 * task: get understand what function of these code? 
 */
function pwy_tool_theme_registry_alter(&$theme_registry) {

  /**
  $theme_registry_copy = $theme_registry;
  $module_path = drupal_get_path('module', 'go_tool') . '/theme';
  _theme_process_registry($theme_registry_copy, 'phptemplate', 'theme_engine', 'my_custom_theme', $module_path);
  $theme_registry += array_diff_key($theme_registry_copy, $theme_registry);
  // A list of templates the module will provide templates for
  $hooks = array('page');
  foreach ($hooks as $hook) {
    // Add the key 'theme paths' if it doesn't exist in this theme's registry
    if (!isset($theme_registry[$hook]['theme paths'])) {
      $theme_registry[$hook]['theme paths'] = array();
    }
    // Shift this module's directory to the top of the theme path list
    if (is_array($theme_registry[$hook]['theme paths'])) {
      $first_element = array_shift($theme_registry[$hook]['theme paths']);
      if ($first_element) {
        array_unshift($theme_registry[$hook]['theme paths'], $first_element, $module_path);
      }
      else {
        array_unshift($theme_registry[$hook]['theme paths'], $module_path);
      }
    }
  }
  */
}

/**
 * Makes the tripal job_id unrecognizable.
 *
 * @param $job_id
 *   The tripal job_id of the analysis you want to make secret.
 *
 * @return
 *   A short string representing the job_id.
 */
function pwy_tool_make_secret($job_id) {

  $mapping = pwy_tool_secret_mapping();
  $secret = str_replace(array_keys($mapping), $mapping, $job_id);

  return $secret;
}

/**
 * Reveals the true job_id for your secret analysis result.
 *
 * @param $secret
 *    The job_id previously made secret by go_tool_make_secret().
 *
 * @return
 *    The revealed tripal job_id.
 */
function pwy_tool_reveal_secret($secret) {

  $mapping = pwy_tool_secret_mapping(TRUE);
  $job_id = str_replace(array_keys($mapping), $mapping, $secret);

  // Check that the job_id exists if it is an integer.
  if (is_numeric($job_id)) {

    $exists = db_query('SELECT job_id FROM {tripal_jobs} WHERE job_id=:id', array(':id' => $job_id))->fetchField();
    if (!$exists) {
      tripal_report_error(
        'pwy_tool',
        TRIPAL_ERROR,
        'Unable to decode the job_id from :id.',
        array(':id' => $secret)
      );
    }
    else {
      return $job_id;
    }
  }
  // Last ditch effort: maybe this job was encoded before the upgrade?
  else {

    $job_id = base64_decode($secret);
    if (is_numeric($job_id)) {
      $exists = db_query('SELECT job_id FROM {tripal_jobs} WHERE job_id=:id', array(':id' => $job_id))->fetchField();
      if (!$exists) {
        tripal_report_error(
          'pwy_tool',
          TRIPAL_ERROR,
          'Unable to decode the job_id from :id.',
          array(':id' => $secret)
        );
      }
      else {
        return $job_id;
      }
    }
    else {
      tripal_report_error(
        'pwy_tool',
        TRIPAL_ERROR,
        'Unable to decode the job_id from :id.',
        array(':id' => $secret)
      );
    }
  }

  return FALSE;
}

/**
 * A single location for keeping track of the mapping used in our secrets.
 */
function pwy_tool_secret_mapping($reveal = FALSE) {
  $mapping = array(
    1 => 'P',
    2 => 'sA',
    3 => 'b',
    4 => 'Q',
    5 => 'Hi',
    6 => 'yt',
    7 => 'f',
    8 => 'zE',
    9 => 'Km',
    0 => 'jVo',
  );

  // Since this is an open-source module with all the code publically available,
  // our secret is not very secret... We are ok with this since the liklihood of
  // profiting by stealing random blast results is pretty low. That said, if this bothers
  // you, feel free to implement the following function in a private module to change
  // this mapping to something that cannot easily be looked up on github. ;-).
  // NOTE: Ensure that the mapping you come up with is unique to ensure that the
  // job_id can be consistently revealed or your users might end up unable to find
  // their own blast results...
  if (function_exists('private_make_mapping_ultra_secret')) {
    $mapping = private_make_mapping_ultra_secret($mapping);
  }

  if ($reveal) {
    return array_flip($mapping);
  }
  else {
    return $mapping;
  }
}

/**
 * Tests your secret mapping over a set of random integers
 * to ensure the job_id can be recovered.
 *
 * @param $num_iterations
 *    An integer representing the number of times you wish to test your mapping.
 */
function pwy_tool_test_secret_mapping($num_iterations = 10000) {

  $all_work = TRUE;

  for($i = 0; $i <= $num_iterations; $i++) {
    $job_id = mt_rand(0,100000);

    $secret = pwy_tool_make_secret($job_id);
    $recovered_job_id = pwy_tool_reveal_secret($secret);

    if ($job_id != $recovered_job_id) {
      drupal_set_message("Unable to recover job_id: $job_id; Secret: $secret.",'error');
      $all_work = FALSE;
    }
  }

  if ($all_work) {
    drupal_Set_message("Secret Mapping works over $num_iterations iterations with random integers.");
  }
}

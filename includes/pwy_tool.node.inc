<?php

/**
 * @file
 * Contains all functions for creating the pwy (Pathway File) node type
 */

/**
 * Implements hook_node_info().
 */
function pwy_tool_node_info() {
  return array(
    'pwy' => array(
      'name' => t('Pathway File'),
      'base' => 'pwy',
      'description' => t('Add a pwy (pathway file) for pathway enrichment analysis.'),
    ),
  );
}

/**
 * Implements hook_node_access().
 */
function pwy_node_access($node, $op, $account) {
  $node_type = $node;
  if (is_object($node)) {
    $node_type = $node->type;
  }

  if($node_type == 'pwy') {
    if ($op == 'create') {
      if (!user_access('create pwy File', $account)) {
        return NODE_ACCESS_DENY;
      }
      return NODE_ACCESS_ALLOW;
    }
    if ($op == 'update') {
      if (!user_access('edit pwy File', $account)) {
        return NODE_ACCESS_DENY;
      }
    }
    if ($op == 'delete') {
      if (!user_access('delete pwy File', $account)) {
        return NODE_ACCESS_DENY;
      }
    }
    if ($op == 'view') {
      if (!user_access('access pwy File', $account)) {
        return NODE_ACCESS_DENY;
      }
    }

    return NODE_ACCESS_IGNORE;
  }
}

/**
 * Form constructor for the pwy node
 *
 * @see pwy_insert()
 * @see pwy_update()
 * @see pwy_delete()
 * @see pwy_load()
 */
function pwy_form($node, &$form_state) {
  $form = array();

  $form['#validate'] = array('pwy_form_validate');
  
  $form['core'] = array(
    '#type' => 'fieldset',
    '#title' => 'General'
  );

  $form['core']['db_name']= array(
    '#type' => 'textfield',
    '#title' => t('Human-readable Name for pathway annotation file'),
    '#required' => TRUE,
    '#default_value' => isset($node->db_name) ? $node->db_name : '',
    '#description' => t('For example: watermelon genome (97103)'),
  );

  $form['core']['db_path']= array(
    '#type' => 'textfield',
    '#title' => t('Full path of pathway annotation file'),
    '#description' => t('For example, /var/www/html/sites/all/modules/tripal_pathway/example/wm_97103_kept.pwy'),
    '#required' => TRUE,
    '#default_value' => isset($node->db_path) ? $node->db_path : '',
  );

  $form['core']['db_idlink'] = array(
    '#type' => 'textfield',
    '#title' => t('Link of gene ID'),
    '#description' => t('For example，www.watermelongenome.org/feature/gene/'),
    '#required' => TRUE,
    '#default_value' => isset($node->db_idlink) ? $node->db_idlink : '',
  );

  $form['core']['db_dblink'] = array(
    '#type' => 'textfield',
    '#title' => t('Link of the pathway ID'),
    '#description' => t('For example, www.watermelongenome.org:1555/WATERMELON_GENOME/new-image?object='),
    '#required' => TRUE,
    '#default_value' => isset($node->db_dblink) ? $node->db_dblink : '',
  );

  $form['core']['db_example_id'] = array(
    '#type' => 'textarea',
    '#title' => t('Example of IDs. One ID for each line '),
    '#required' => TRUE,
    '#default_value' => isset($node->db_example_id) ? $node->db_example_id : '',
    '#description' => t('The IDs will be used as example for correspond GAF.'),
    '#cols' => 30,
    '#rows' => 10,
    '#maxlength' => 1000,
  );

  return $form;
}

function pwy_form_validate($form, $form_state) {

}

/**
 * Implements hook_insert().
 */
function pwy_insert($node) {
  
  // Actually insert the record.
  db_insert('pwy')->fields(array(
    'nid'                 => $node->nid,
    'name'                => $node->db_name,
    'path'                => $node->db_path,
    'idlink'              => $node->db_idlink,
    'dblink'              => $node->db_dblink,
    'idexample'           => $node->db_example_id,
  ))->execute();
}

/**
 * Implements hook_node_insert().
 * This function acts on ALL NODES
 */
function pwy_tool_node_insert($node) {
  if ($node->type == 'pwy') {
    $node->title = $node->db_name;
  }
}

/**
 * Implements hook_update().
 * modify pwy in database level
 */
function pwy_update($node) {
  // Update the record.
  db_update('pwy')->fields(array(
    'name'                => $node->db_name,
    'path'                => $node->db_path,
    'idlink'              => $node->db_idlink,
    'dblink'              => $node->db_dblink,
    'idexample'           => $node->db_example_id,
  ))->condition('nid', $node->nid)->execute();
}

/**
 * Implements hook_node_update().
 * This function acts on ALL NODES
 * for pwy_tool_node
 */
function pwy_tool_node_update($node) {
  if ($node->type == 'pwy')  {
    $node->title = $node->db_name;
  }
}

/**
 * Implements hook_delete(). 
 * modify pwy in database level 
 */
function pwy_delete($node) {
  db_delete('pwy')->condition('nid',$node->nid)->execute();
}

/**
 * Implements hook_load().
 * modify pwy in database level
 */
function pwy_load($nodes) {

  $sql = "
    SELECT nid, name, path, idlink, dblink, idexample
    FROM {pwy} 
    WHERE nid IN (:nids)";
  $result = db_query($sql, array(':nids' => array_keys($nodes)));

  foreach ($result as $record) {
    $nodes[$record->nid]->db_name = $record->name;
    $nodes[$record->nid]->db_path = $record->path;
    $nodes[$record->nid]->title = $record->name;
    $nodes[$record->nid]->db_idlink = $record->idlink;
    $nodes[$record->nid]->db_dblink = $record->dblink;
    $nodes[$record->nid]->db_example_id = $record->idexample;
  }
}


# DNA

Version：

<table>
  <tr>
    <th rowspan="2">Model</th>
    <th rowspan="2">Description</th>
    <th rowspan="2">Final feature processing</th>
    <th colspan="4">Parameter</th>
    <th rowspan="2">Threshold</th>
    <th rowspan="2">FOLD</th>
  </tr>

  <tr>
    <th>Normalization</th>
    <th>Dropout</th>
    <th>Learning Rate</th>
    <th>Weight</th>
  </tr>

  <tr>
    <td >V1</td>
    <td> Baseline </td>
    <td rowspan="3">MaxPool → Global Average Pooling → Flatten</td>
    <td>BatchNorm</td>
    <td>0</td>
    <td>0.0001</td>
    <td>0.0001</td>
    <td rowspan="6">0.5</td>
    <td rowspan="6">5-fold</td>
  </tr>

  <tr>
    <td >V2</td>
    <td> Stronger regularization </td>
    <td>BatchNorm</td>
    <td rowspan="5">0.5</td>
    <td rowspan="5">0.00005</td>
    <td rowspan="5">0.001</td>
  </tr>

  <tr>
    <td >V3</td>
    <td> Replaced BatchNorm with GroupNorm </td>
    <td rowspan="4">GroupNorm</td>
  </tr>

  <tr>
    <td >V4</td>
    <td> Replaced GAP with Average Pooling  </td>
    <td >MaxPool → AvgPool → Flatten</td>
  </tr>

  <tr>
    <td >V5</td>
    <td> Remove AvgPool  </td>
    <td >MaxPool → Flatten</td>
  </tr>

  <tr>
    <td >V6</td>
    <td> Remove MaxPool   </td>
    <td >Global Average Pooling → Flatten</td>
  </tr>
  
<table>
 

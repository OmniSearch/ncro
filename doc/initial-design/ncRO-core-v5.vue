<!-- Tufts VUE 3.2.2 concept-map (ncRO-core-v5.vue) 2015-06-08 -->
<!-- Tufts VUE: http://vue.tufts.edu/ -->
<!-- Do Not Remove: VUE mapping @version(1.1) jar:file:/C:/Program%20Files/VUE/VUE.jar!/tufts/vue/resources/lw_mapping_1_1.xml -->
<!-- Do Not Remove: Saved date Mon Jun 08 23:18:28 CDT 2015 by Jingshan Huang on platform Windows 7 6.1 in JVM 1.7.0_21-b11 -->
<!-- Do Not Remove: Saving version @(#)VUE: built May 23 2013 at 2146 by tomadm on Linux 2.6.18-348.2.1.el5 i386 JVM 1.7.0_21-b11(bits=32) -->
<?xml version="1.0" encoding="US-ASCII"?>
<LW-MAP xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:noNamespaceSchemaLocation="none" ID="0" label="ncRO-core-v5.vue"
    created="1431734769931" x="0.0" y="0.0" width="1.4E-45"
    height="1.4E-45" strokeWidth="0.0" autoSized="false">
    <resource referenceCreated="1433823508676" size="59718"
        spec="C:\Users\Jingshan Huang\Documents\GitHub\OmniSearch\ncRO-Ontology-Design\ncRO-core-v5.vue"
        type="1" xsi:type="URLResource">
        <title>ncRO-core-v5.vue</title>
        <property key="File" value="C:\Users\Jingshan Huang\Documents\GitHub\OmniSearch\ncRO-Ontology-Design\ncRO-core-v5.vue"/>
    </resource>
    <fillColor>#FFFFFF</fillColor>
    <strokeColor>#404040</strokeColor>
    <textColor>#000000</textColor>
    <font>SansSerif-plain-14</font>
    <URIString>http://vue.tufts.edu/rdf/resource/5a6f927f0a147a5701bb3ae4d29d2d28</URIString>
    <child ID="26" label="SO: ncRNA" layerID="1" created="1431736225350"
        x="312.4991" y="532.7976" width="72.0" height="23.0"
        strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#9DDB53</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/5a6f92980a147a5701bb3ae40b70afce</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="28" label="SO: small_regulatory_ncRNA" layerID="1"
        created="1431736262255" x="432.72015" y="489.4687" width="171.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#9DDB53</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/5a6f929a0a147a5701bb3ae44053001c</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="29" label="is_a" layerID="1" created="1431736282909"
        x="381.14386" y="511.9687" width="92.53024" height="24.367188"
        strokeWidth="1.0" autoSized="false" controlCount="0"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/5a6f929b0a147a5701bb3ae45b913737</URIString>
        <point1 x="473.17413" y="512.4687"/>
        <point2 x="381.6439" y="535.8359"/>
        <ID1 xsi:type="node">28</ID1>
        <ID2 xsi:type="node">26</ID2>
    </child>
    <child ID="30" label="SO: miRNA" layerID="1" created="1431736299454"
        x="721.4205" y="324.5058" width="73.0" height="23.0"
        strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#9DDB53</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/5a6f929d0a147a5701bb3ae4fb5fe718</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="31" label="is_a" layerID="1" created="1431736316435"
        x="527.91815" y="347.0058" width="205.12164" height="142.96289"
        strokeWidth="1.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/5a6f929f0a147a5701bb3ae4d648178b</URIString>
        <point1 x="732.5398" y="347.5058"/>
        <point2 x="528.41815" y="489.4687"/>
        <ID1 xsi:type="node">30</ID1>
        <ID2 xsi:type="node">28</ID2>
        <ctrlPoint0 x="601.76105" y="406.76163" xsi:type="point"/>
    </child>
    <child ID="46" label="SO: gene" layerID="1" created="1431740614381"
        x="973.59924" y="326.1799" width="62.0" height="23.0"
        strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#9DDB53</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b00a147a5701bb3ae44cf956bd</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="47" label="ncro: miRNA_target_gene" layerID="1"
        created="1431740720414" x="620.79315" y="523.52997"
        width="153.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b00a147a5701bb3ae467b66d03</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="48" label="is_a" layerID="1" created="1431740743700"
        x="772.9625" y="348.67993" width="227.18457" height="189.98761"
        strokeWidth="1.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b10a147a5701bb3ae465962ea4</URIString>
        <point1 x="773.4625" y="538.16754"/>
        <point2 x="999.6471" y="349.17993"/>
        <ID1 xsi:type="node">47</ID1>
        <ID2 xsi:type="node">46</ID2>
        <ctrlPoint0 x="915.74084" y="544.02826" xsi:type="point"/>
    </child>
    <child ID="59" label="ncro: miRNA_and_target_gene_binding"
        layerID="1" created="1433381847353" x="627.717" y="426.31097"
        width="229.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc3bc18dc0a80c610191f5df11a59903</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="60" label="ncro: is_about_miRNA_target_gene" layerID="1"
        created="1433381966836" x="629.9612" y="448.81104" width="171.0"
        height="73.15405" strokeWidth="1.0" autoSized="false"
        controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#7C18C9</strokeColor>
        <textColor>#7C18C9</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc3bc18ec0a80c610191f5dfe25d1dd6</URIString>
        <point1 x="733.0222" y="449.31104"/>
        <point2 x="715.4742" y="521.4651"/>
        <ID1 xsi:type="node">59</ID1>
        <ctrlPoint0 x="706.6742" y="482.26508" xsi:type="point"/>
    </child>
    <child ID="61" label="ncro: is_about_miRNA" layerID="1"
        created="1433382183779" x="698.0063" y="348.09583" width="107.0"
        height="78.71518" strokeWidth="1.0" autoSized="false"
        controlCount="0" arrowState="2" xsi:type="link">
        <strokeColor>#7C18C9</strokeColor>
        <textColor>#7C18C9</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc3cec2dc0a80c610191f5df3bd3b19f</URIString>
        <point1 x="744.3384" y="426.31104"/>
        <point2 x="758.6742" y="348.59586"/>
        <ID1 xsi:type="node">59</ID1>
    </child>
    <child ID="63" label="ncro: human_miRNA" layerID="1"
        created="1433382370738" x="61.92553" y="305.7788" width="126.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc40e990c0a80c610191f5df36ddb9b2</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="64" label="ncro: mouse_miRNA" layerID="1"
        created="1433382395274" x="85.35077" y="355.12073" width="126.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc40e990c0a80c610191f5dfce5c1ff5</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="65" label="is_a" layerID="1" created="1433382425225"
        x="210.33228" y="366.48682" width="86.94035" height="14.0"
        strokeWidth="1.0" autoSized="false" controlCount="0"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc40e990c0a80c610191f5df72742021</URIString>
        <point1 x="210.83229" y="370.68896"/>
        <point2 x="296.77264" y="376.28467"/>
        <ID1 xsi:type="node">64</ID1>
        <ID2 xsi:type="node">91</ID2>
    </child>
    <child ID="66" label="is_a" layerID="1" created="1433382437823"
        x="180.7211" y="327.61346" width="187.6109" height="42.63919"
        strokeWidth="1.0" autoSized="false" controlCount="2"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc40e991c0a80c610191f5dfa3a027d6</URIString>
        <point1 x="181.2211" y="328.11346"/>
        <point2 x="367.832" y="369.75266"/>
        <ID1 xsi:type="node">63</ID1>
        <ID2 xsi:type="node">91</ID2>
        <ctrlPoint0 x="275.47247" y="346.2531" xsi:type="point"/>
        <ctrlPoint1 x="346.69827" y="323.3757" xsi:type="point"/>
    </child>
    <child ID="68" layerID="1" created="1433382524167" x="-106.42582"
        y="256.13504" width="150.75" height="162.25" strokeWidth="1.0"
        autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc4a5781c0a80c610191f5df737bb488</URIString>
        <child ID="72" label="ncro: hsa-miR-125b"
            created="1433382574535" x="34.0" y="22.0" width="148.0"
            height="27.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#FDE888</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/bc4a5782c0a80c610191f5df230e2262</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <child ID="74" label="ncro: hsa-let-7a-1"
            created="1433382849934" x="34.0" y="45.25" width="132.0"
            height="27.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#FDE888</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/bc4a5782c0a80c610191f5dfc411e3c4</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <child ID="76" label="ncro: hsa-miR-10b" created="1433382974053"
            x="34.0" y="68.5" width="139.0" height="27.0"
            strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#FDE888</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/bc4a5782c0a80c610191f5dfd082b240</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <child ID="77" label="ncro: hsa-miR-21" created="1433383060372"
            x="34.0" y="91.75" width="130.0" height="27.0"
            strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#FDE888</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/bc4a5783c0a80c610191f5dfcec5ba4f</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <child ID="81" label="ncro: hsa-miR-30e" created="1433383413170"
            x="34.0" y="115.0" width="139.0" height="27.0"
            strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#FDE888</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/bc519b59c0a80c6101fc8aab6accbf77</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <child ID="82" label="And so on..." created="1433383440561"
            x="34.0" y="138.25" width="102.0" height="18.0"
            strokeWidth="0.0" autoSized="false" xsi:type="text">
            <strokeColor>#404040</strokeColor>
            <textColor>#000000</textColor>
            <font>-plain-1</font>
            <URIString>http://vue.tufts.edu/rdf/resource/bc519b5ac0a80c6101fc8aabed32d421</URIString>
            <richText>&lt;html&gt;
  &lt;head style="color: #000000" color="#000000"&gt;
    &lt;style type="text/css"&gt;
      &lt;!--
        body { font-family: Arial; margin-bottom: 0px; margin-right: 0px; font-size: 12; margin-top: 0px; margin-left: 0px; color: #000000 }
        ol { font-family: Arial; font-size: 12; margin-top: 6; margin-left: 30; list-style-position: outside; vertical-align: middle }
        p { margin-bottom: 0; margin-right: 0; margin-top: 0; margin-left: 0; color: #000000 }
        ul { font-family: Arial; font-size: 12; margin-top: 6; margin-left: 30; list-style-position: outside; vertical-align: middle }
      --&gt;
    &lt;/style&gt;
    
  &lt;/head&gt;
  &lt;body&gt;
    &lt;p style="color: #000000" color="#000000"&gt;
      And so on...
    &lt;/p&gt;
  &lt;/body&gt;
&lt;/html&gt;
</richText>
            <label>And so on...</label>
        </child>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="75" label="is_a" layerID="1" created="1433382954228"
        x="43.82418" y="277.76074" width="74.73177" height="28.518066"
        strokeWidth="1.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc4a5783c0a80c610191f5df2f213378</URIString>
        <point1 x="44.32418" y="296.36063"/>
        <point2 x="118.055954" y="305.7788"/>
        <ID1 xsi:type="node">68</ID1>
        <ID2 xsi:type="node">63</ID2>
        <ctrlPoint0 x="95.75846" y="268.45175" xsi:type="point"/>
    </child>
    <child ID="88" label="is_a" layerID="1" created="1433548171824"
        x="767.3627" y="246.09421" width="129.87634" height="180.71677"
        strokeWidth="1.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c621d82ec0a80c4a00148dd87f0103cd</URIString>
        <point1 x="767.8627" y="426.31097"/>
        <point2 x="856.37134" y="246.59421"/>
        <ID1 xsi:type="node">59</ID1>
        <ID2 xsi:type="node">132</ID2>
        <ctrlPoint0 x="962.3611" y="339.0942" xsi:type="point"/>
    </child>
    <child ID="91" label="ncro: mammalian_miRNA" layerID="1"
        created="1433548846994" x="296.07248" y="369.75262"
        width="154.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c63bb3b9c0a80c4a00148dd84aae1949</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="92" label="ncro: non_mammalian_miRNA" layerID="1"
        created="1433548983258" x="288.87247" y="433.7526" width="182.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c63bb3b9c0a80c4a00148dd88dc3cf22</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="100" label="is_a" layerID="1" created="1433549793675"
        x="408.64722" y="335.10516" width="313.27325" height="35.14746"
        strokeWidth="1.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c63bb3bdc0a80c4a00148dd88ff8e05c</URIString>
        <point1 x="409.14722" y="369.75262"/>
        <point2 x="721.4205" y="335.60516"/>
        <ID1 xsi:type="node">91</ID1>
        <ID2 xsi:type="node">30</ID2>
        <ctrlPoint0 x="523.0944" y="333.42828" xsi:type="point"/>
    </child>
    <child ID="101" label="is_a" layerID="1" created="1433549848124"
        x="411.2212" y="342.80414" width="312.77014" height="91.448456"
        strokeWidth="1.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c63bb3bec0a80c4a00148dd827b3c8ab</URIString>
        <point1 x="411.7212" y="433.7526"/>
        <point2 x="723.49133" y="343.30414"/>
        <ID1 xsi:type="node">92</ID1>
        <ID2 xsi:type="node">30</ID2>
        <ctrlPoint0 x="575.0944" y="374.76163" xsi:type="point"/>
    </child>
    <child ID="102" label="ncro: miRNA_expression" layerID="1"
        created="1433550036481" x="292.75687" y="299.1521" width="148.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c63ed923c0a80c4a00148dd824409ce8</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="103" label="RO: has_participant" layerID="1"
        created="1433550220480" x="400.53485" y="288.04312"
        width="323.37347" height="41.287323" strokeWidth="1.0"
        autoSized="false" controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c640fac5c0a80c4a00148dd82c12f016</URIString>
        <point1 x="401.03485" y="299.1521"/>
        <point2 x="723.4083" y="328.83044"/>
        <ID1 xsi:type="node">102</ID1>
        <ID2 xsi:type="node">30</ID2>
        <ctrlPoint0 x="469.76108" y="276.09494" xsi:type="point"/>
    </child>
    <child ID="104" label="GO: transcription" layerID="1"
        created="1433550288241" x="868.6929" y="146.95157" width="103.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#B5B995</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c644475ac0a80c4a00148dd894a8124b</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="105" label="GO: translation" layerID="1"
        created="1433550311808" x="245.09164" y="192.01862" width="93.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#B5B995</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c644475ac0a80c4a00148dd8c44e24cb</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="106" label="UBERON: tissue" layerID="1"
        created="1433550320530" x="70.29166" y="192.15195" width="103.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#F2AE45</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c644475ac0a80c4a00148dd8ef2fd2a0</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="107" label="PRO: protein" layerID="1"
        created="1433550323639" x="107.582016" y="75.28398" width="81.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#83CEFF</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c644475bc0a80c4a00148dd8352bab0c</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="108" label="ncro: miRNA_expressed_in_tissue" layerID="1"
        created="1433550673037" x="59.435852" y="214.65195"
        width="239.13939" height="86.318085" strokeWidth="1.0"
        autoSized="false" controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#7C18C9</strokeColor>
        <textColor>#7C18C9</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c647d40dc0a80c4a00148dd83036bf32</URIString>
        <point1 x="298.07523" y="300.47003"/>
        <point2 x="114.27931" y="215.15195"/>
        <ID1 xsi:type="node">102</ID1>
        <ID2 xsi:type="node">106</ID2>
        <ctrlPoint0 x="79.69443" y="268.09494" xsi:type="point"/>
    </child>
    <child ID="111" label="ncro: regulate_mRNA_translation" layerID="1"
        created="1433550883004" x="257.0208" y="214.51855" width="159.0"
        height="85.133545" strokeWidth="1.0" autoSized="false"
        controlCount="2" arrowState="2" xsi:type="link">
        <strokeColor>#7C18C9</strokeColor>
        <textColor>#7C18C9</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c64b75a4c0a80c4a00148dd84ea32033</URIString>
        <point1 x="369.40402" y="299.1521"/>
        <point2 x="292.63034" y="215.01855"/>
        <ID1 xsi:type="node">102</ID1>
        <ID2 xsi:type="node">105</ID2>
        <ctrlPoint0 x="379.79114" y="254.02724" xsi:type="point"/>
        <ctrlPoint1 x="296.91953" y="262.50452" xsi:type="point"/>
    </child>
    <child ID="112" label="RO: has_output" layerID="1"
        created="1433715059634" x="139.5052" y="98.49585"
        width="127.76497" height="94.022766" strokeWidth="1.0"
        autoSized="false" controlCount="2" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d0158ec6c0a80c4a019878a2bb93f77a</URIString>
        <point1 x="266.77017" y="192.01862"/>
        <point2 x="169.07419" y="98.99585"/>
        <ID1 xsi:type="node">105</ID1>
        <ctrlPoint0 x="147.2484" y="136.64314" xsi:type="point"/>
        <ctrlPoint1 x="179.48405" y="152.346" xsi:type="point"/>
    </child>
    <child ID="113" layerID="1" created="1433715286761" x="-97.76742"
        y="-108.963165" width="183.0" height="109.0" strokeWidth="1.0"
        autoSized="true" xsi:type="node">
        <fillColor>#00000000</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d019f1f3c0a80c4a019878a2cb4655b9</URIString>
        <child ID="114" label="UBERON: organ" created="1433715286761"
            x="34.0" y="22.0" width="126.0" height="27.0"
            strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#F2AE45</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/d019f1f3c0a80c4a019878a2ad6ee2a1</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <child ID="115" label="NCIt: Organ" created="1433715286761"
            x="34.0" y="45.25" width="92.0" height="27.0"
            strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#DAA9FF</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/d019f1f3c0a80c4a019878a277b1059b</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <child ID="116"
            label="SNOMEDCT: Anatomical&#xa;or acquired body structure"
            created="1433715286761" x="34.0" y="68.5" width="191.0"
            height="46.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#8AEE95</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/d019f1f3c0a80c4a019878a2012220a0</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="120" label="BFO2: occurs_in" layerID="1"
        created="1433715499187" x="67.96263" y="-6.3041506" width="81.0"
        height="82.088135" strokeWidth="1.0" autoSized="false"
        controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d01ba71fc0a80c4a019878a297eb5ded</URIString>
        <point1 x="139.12088" y="75.28398"/>
        <point2 x="85.87418" y="-5.8041506"/>
        <ID1 xsi:type="node">107</ID1>
        <ctrlPoint0 x="104.427734" y="30.761618" xsi:type="point"/>
    </child>
    <child ID="121" layerID="1" created="1433715637865" x="155.97443"
        y="-152.73811" width="173.25" height="71.5" strokeWidth="1.0"
        autoSized="true" xsi:type="node">
        <fillColor>#00000000</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d01d6b7fc0a80c4a019878a25c28fe32</URIString>
        <child ID="123" label="NCIt: Organism" created="1433715637875"
            x="34.0" y="22.0" width="117.0" height="27.0"
            strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#DAA9FF</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/d01d6b7fc0a80c4a019878a2207b15e2</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <child ID="124" label="SNOMEDCT: Organism"
            created="1433715637875" x="34.0" y="45.25" width="178.0"
            height="27.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#8AEE95</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/d01d6b7fc0a80c4a019878a2850a5cf0</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="125" layerID="1" created="1433715751553" x="403.40747"
        y="-182.07172" width="183.0" height="101.5" strokeWidth="1.0"
        autoSized="false" xsi:type="node">
        <fillColor>#00000000</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d01f8535c0a80c4a019878a2e361ecd0</URIString>
        <child ID="126" label="DOID: disease" created="1433715751553"
            x="34.0" y="22.0" width="114.0" height="27.0"
            strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#C1F780</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/d01f8535c0a80c4a019878a268cc4552</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <child ID="127" label="NCIt: Disease or Disorder"
            created="1433715751553" x="34.0" y="45.25" width="191.0"
            height="27.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#DAA9FF</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/d01f8535c0a80c4a019878a2fb671f54</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <child ID="128" label="SNOMEDCT: Disease"
            created="1433715751553" x="34.0" y="68.5" width="169.0"
            height="27.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#8AEE95</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/d01f8535c0a80c4a019878a29242d437</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="129" label="RO: part_of" layerID="1"
        created="1433715865314" x="-6.3017015" y="-154.88974"
        width="159.87589" height="46.426605" strokeWidth="1.0"
        autoSized="false" controlCount="2" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d020a786c0a80c4a019878a20d20e95a</URIString>
        <point1 x="-5.8017015" y="-108.963135"/>
        <point2 x="153.07419" y="-140.20415"/>
        <ID1 xsi:type="node">113</ID1>
        <ctrlPoint0 x="-5.3445253" y="-162.46358" xsi:type="point"/>
        <ctrlPoint1 x="85.28053" y="-148.8533" xsi:type="point"/>
    </child>
    <child ID="130" label="BFO2: has_disposition" layerID="1"
        created="1433715928785" x="214.22266" y="-195.40414"
        width="188.15152" height="43.166046" strokeWidth="1.0"
        autoSized="false" controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d021a598c0a80c4a019878a2a4413832</URIString>
        <point1 x="234.20139" y="-152.7381"/>
        <point2 x="401.87418" y="-160.9349"/>
        <ID1 xsi:type="node">121</ID1>
        <ctrlPoint0 x="218.4075" y="-219.97177" xsi:type="point"/>
    </child>
    <child ID="131" label="RO: participates_in" layerID="1"
        created="1433750311138" x="925.412" y="169.45154" width="91.0"
        height="157.2284" strokeWidth="1.0" autoSized="false"
        controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2325cc2c0a80c4a00cfd279cd814c0c</URIString>
        <point1 x="1001.36115" y="326.17993"/>
        <point2 x="928.23126" y="169.95154"/>
        <ID1 xsi:type="node">46</ID1>
        <ID2 xsi:type="node">104</ID2>
        <ctrlPoint0 x="977.0277" y="239.76085" xsi:type="point"/>
    </child>
    <child ID="132" label="GO: binding" layerID="1"
        created="1433750793012" x="805.1942" y="223.5942" width="76.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#B5B995</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2356455c0a80c4a00cfd27986d5571f</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="135" label="CLO: cell line" layerID="1"
        created="1433751135911" x="942.3611" y="-74.40573" width="84.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#BDE5F2</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d23e9359c0a80c4a00cfd2791b229471</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="137" label="ncro: animal_study_model" layerID="1"
        created="1433751259326" x="833.6932" y="40.927067" width="157.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d23e9359c0a80c4a00cfd2797af32b2d</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="138" label="is_a" layerID="1" created="1433751404809"
        x="329.37418" y="-89.50415" width="550.52246" height="130.93121"
        strokeWidth="1.0" autoSized="false" controlCount="2"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d23e9359c0a80c4a00cfd2797ac1d191</URIString>
        <point1 x="879.39667" y="40.927063"/>
        <point2 x="329.87418" y="-89.00415"/>
        <ID1 xsi:type="node">137</ID1>
        <ctrlPoint0 x="409.7187" y="-123.764145" xsi:type="point"/>
        <ctrlPoint1 x="615.4363" y="0.07449722" xsi:type="point"/>
    </child>
    <child ID="139" label="BFO2: occurs_in" layerID="1"
        created="1433751442012" x="879.02124" y="64.165085" width="81.0"
        height="83.28804" strokeWidth="1.0" autoSized="false"
        controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d23f688cc0a80c4a00cfd279de4fe6ee</URIString>
        <point1 x="920.28845" y="146.95312"/>
        <point2 x="916.27423" y="64.665085"/>
        <ID1 xsi:type="node">104</ID1>
        <ctrlPoint0 x="920.76105" y="90.094955" xsi:type="point"/>
    </child>
    <child ID="140" label="BFO2: occurs_in" layerID="1"
        created="1433751458806" x="931.1273" y="-51.034912"
        width="142.16223" height="198.48645" strokeWidth="1.0"
        autoSized="false" controlCount="2" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d23f688cc0a80c4a00cfd2790718175d</URIString>
        <point1 x="931.6273" y="146.95154"/>
        <point2 x="1004.27423" y="-50.534912"/>
        <ID1 xsi:type="node">104</ID1>
        <ctrlPoint0 x="1080.5315" y="-2.8061905" xsi:type="point"/>
        <ctrlPoint1 x="1028.2734" y="-5.7365265" xsi:type="point"/>
    </child>
    <child ID="142" label="ncro: is_model_of_disease" layerID="1"
        created="1433751602281" x="586.9742" y="-111.834915"
        width="311.8158" height="153.262" strokeWidth="1.0"
        autoSized="false" controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#7C18C9</strokeColor>
        <textColor>#7C18C9</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d24213edc0a80c4a00cfd2792c02a4fd</URIString>
        <point1 x="898.29" y="40.927063"/>
        <point2 x="587.4742" y="-111.334915"/>
        <ID1 xsi:type="node">137</ID1>
        <ctrlPoint0 x="755.02765" y="-77.57248" xsi:type="point"/>
    </child>
    <child ID="143" label="ncro: is_model_of_disease" layerID="1"
        created="1433751632123" x="586.9742" y="-135.03491"
        width="363.32123" height="61.647125" strokeWidth="1.0"
        autoSized="false" controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#7C18C9</strokeColor>
        <textColor>#7C18C9</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d24213edc0a80c4a00cfd279c290fa2e</URIString>
        <point1 x="949.7954" y="-73.88779"/>
        <point2 x="587.4742" y="-134.53491"/>
        <ID1 xsi:type="node">135</ID1>
        <ctrlPoint0 x="753.87415" y="-136.13492" xsi:type="point"/>
    </child>
    <child ID="144" label="ncro: transcription_of_miRNA" layerID="1"
        created="1433751984266" x="447.6946" y="16.260838" width="172.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d24bfccfc0a80c4a00cfd27909ff6308</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="145" label="is_a" layerID="1" created="1433752005899"
        x="559.8956" y="-2.2260985" width="336.0786" height="148.1912"
        strokeWidth="1.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d24bfccfc0a80c4a00cfd27919c338fa</URIString>
        <point1 x="560.3956" y="16.260849"/>
        <point2 x="895.47424" y="145.46509"/>
        <ID1 xsi:type="node">144</ID1>
        <ctrlPoint0 x="721.76105" y="-53.238384" xsi:type="point"/>
    </child>
    <child ID="147" label="ncro: miRNA_transcription_initiation"
        layerID="1" created="1433752577010" x="354.361" y="103.59421"
        width="208.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d25561dac0a80c4a00cfd27933737d8f</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="148" label="ncro: protein_miRNA_promoter_binding"
        layerID="1" created="1433752594051" x="451.6944" y="193.59422"
        width="230.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d25561dac0a80c4a00cfd279ddf26489</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="149" label="ncro: promoter_of_miRNA" layerID="1"
        created="1433752659308" x="645.02765" y="102.26087"
        width="154.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d25561dac0a80c4a00cfd279322bbbff</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="150" label="is_a" layerID="1" created="1433753078241"
        x="586.81244" y="216.09422" width="221.4054" height="41.20117"
        strokeWidth="1.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d25de720c0a80c4a00cfd27961f54cc2</URIString>
        <point1 x="587.31244" y="216.59422"/>
        <point2 x="807.71783" y="243.0641"/>
        <ID1 xsi:type="node">148</ID1>
        <ID2 xsi:type="node">132</ID2>
        <ctrlPoint0 x="684.42773" y="270.76163" xsi:type="point"/>
    </child>
    <child ID="151" label="RO: participates_in" layerID="1"
        created="1433753378695" x="184.07683" y="95.58069"
        width="265.4973" height="107.91516" strokeWidth="1.0"
        autoSized="false" controlCount="2" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d25de720c0a80c4a00cfd27921d628da</URIString>
        <point1 x="184.57683" y="96.08069"/>
        <point2 x="449.07416" y="202.99585"/>
        <ID1 xsi:type="node">107</ID1>
        <ctrlPoint0 x="434.23438" y="159.67857" xsi:type="point"/>
        <ctrlPoint1 x="349.56738" y="195.20567" xsi:type="point"/>
    </child>
    <child ID="152" label="ncro: regulate_miRNA_transcription"
        layerID="1" created="1433753551626" x="151.07207" y="-45.447575"
        width="377.7021" height="121.231575" strokeWidth="1.0"
        autoSized="false" controlCount="2" arrowState="2" xsi:type="link">
        <strokeColor>#AF55F4</strokeColor>
        <textColor>#AF55F4</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2615c64c0a80c4a00cfd27950d3405c</URIString>
        <point1 x="151.57207" y="75.284"/>
        <point2 x="528.2742" y="14.99585"/>
        <ID1 xsi:type="node">107</ID1>
        <ctrlPoint0 x="201.49495" y="-89.215996" xsi:type="point"/>
        <ctrlPoint1 x="487.87418" y="-43.40415" xsi:type="point"/>
    </child>
    <child ID="154" label="RO: precedes" layerID="1"
        created="1433753916773" x="444.17715" y="127.29585"
        width="83.408966" height="66.798386" strokeWidth="1.0"
        autoSized="false" controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d26514fcc0a80c4a00cfd279d6be07aa</URIString>
        <point1 x="527.0861" y="193.59424"/>
        <point2 x="465.07416" y="127.79585"/>
        <ID1 xsi:type="node">148</ID1>
        <ctrlPoint0 x="460.27417" y="174.19585" xsi:type="point"/>
    </child>
    <child ID="155" label="RO: precedes" layerID="1"
        created="1433753966436" x="461.25687" y="40.095848"
        width="79.590515" height="63.99839" strokeWidth="1.0"
        autoSized="false" controlCount="2" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2659098c0a80c4a00cfd279740792b6</URIString>
        <point1 x="461.75687" y="103.59424"/>
        <point2 x="534.6742" y="40.595848"/>
        <ID1 xsi:type="node">147</ID1>
        <ctrlPoint0 x="479.77716" y="42.568653" xsi:type="point"/>
        <ctrlPoint1 x="539.6721" y="109.882" xsi:type="point"/>
    </child>
    <child ID="156" label="RO: participates_in" layerID="1"
        created="1433754046535" x="621.37415" y="27.295849"
        width="148.93018" height="75.46501" strokeWidth="1.0"
        autoSized="false" controlCount="2" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2693664c0a80c4a00cfd279f89efdba</URIString>
        <point1 x="726.40967" y="102.260864"/>
        <point2 x="621.87415" y="27.795849"/>
        <ID1 xsi:type="node">149</ID1>
        <ctrlPoint0 x="742.086" y="61.121044" xsi:type="point"/>
        <ctrlPoint1 x="741.29755" y="27.877655" xsi:type="point"/>
    </child>
    <child ID="157" label="RO: participates_in" layerID="1"
        created="1433754091599" x="550.9742" y="57.89557"
        width="151.54376" height="44.865265" strokeWidth="1.0"
        autoSized="false" controlCount="2" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2693664c0a80c4a00cfd2795daa83d2</URIString>
        <point1 x="702.01794" y="102.26083"/>
        <point2 x="551.4742" y="102.19585"/>
        <ID1 xsi:type="node">149</ID1>
        <ctrlPoint0 x="636.5032" y="64.60807" xsi:type="point"/>
        <ctrlPoint1 x="564.17206" y="40.294544" xsi:type="point"/>
    </child>
    <child ID="158" label="GO: RNA_metabolic_process" layerID="1"
        created="1433754139438" x="229.19443" y="24.594929"
        width="173.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#B5B995</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2693664c0a80c4a00cfd2792b44e1da</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="159" label="is_a" layerID="1" created="1433754164170"
        x="323.77417" y="48.095848" width="74.825134" height="55.99836"
        strokeWidth="1.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2693664c0a80c4a00cfd2793c92f5c8</URIString>
        <point1 x="398.0993" y="103.59421"/>
        <point2 x="324.27417" y="48.595848"/>
        <ID1 xsi:type="node">147</ID1>
        <ctrlPoint0 x="336.27417" y="91.79585" xsi:type="point"/>
    </child>
    <child ID="160" label="SO: promoter" layerID="1"
        created="1433754251067" x="713.1944" y="185.26158" width="84.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#9DDB53</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d26a4329c0a80c4a00cfd2790005f5d0</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="161" label="is_a" layerID="1" created="1433754268989"
        x="728.72845" y="124.760864" width="24.333862" height="61.00061"
        strokeWidth="1.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d26a4329c0a80c4a00cfd279d125459a</URIString>
        <point1 x="729.22845" y="125.260864"/>
        <point2 x="752.5623" y="185.26147"/>
        <ID1 xsi:type="node">149</ID1>
        <ID2 xsi:type="node">160</ID2>
        <ctrlPoint0 x="744.3611" y="149.42825" xsi:type="point"/>
    </child>
    <child ID="162" label="GO: metabolic_process" layerID="1"
        created="1433756364341" x="487.07422" y="240.98132"
        width="141.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#B5B995</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d28adaccc0a80c4a017b3395e1c17d6c</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="163" label="is_a" layerID="1" created="1433756389657"
        x="382.96927" y="255.63687" width="105.04016" height="44.01523"
        strokeWidth="1.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d28adacdc0a80c4a017b339557150a81</URIString>
        <point1 x="383.46927" y="299.1521"/>
        <point2 x="487.50943" y="256.13687"/>
        <ID1 xsi:type="node">102</ID1>
        <ID2 xsi:type="node">162</ID2>
        <ctrlPoint0 x="442.5742" y="258.48132" xsi:type="point"/>
    </child>
    <child ID="164" label="RO: participates_in" layerID="1"
        created="1433762567433" x="561.47095" y="123.25076" width="91.0"
        height="69.84509" strokeWidth="1.0" autoSized="false"
        controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2e8dd3ac0a80c4a002b8aa654caad21</URIString>
        <point1 x="650.061" y="123.75076"/>
        <point2 x="574.6742" y="192.59586"/>
        <ID1 xsi:type="node">149</ID1>
        <ctrlPoint0 x="601.5742" y="130.48134" xsi:type="point"/>
    </child>
    <child ID="172" layerID="1" created="1433778917213" x="61.204132"
        y="413.35632" width="174.75" height="139.0" strokeWidth="1.0"
        autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d3e2e702c0a80c4a00de5ab4b0742fba</URIString>
        <child ID="173" label="ncro: fly_miRNA" created="1433778917213"
            x="34.0" y="22.0" width="122.0" height="27.0"
            strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#FDE888</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/d3e2e703c0a80c4a00de5ab4b1c84610</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <child ID="174" label="ncro: c_elegans_miRNA"
            created="1433778917213" x="34.0" y="45.25" width="180.0"
            height="27.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#FDE888</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/d3e2e704c0a80c4a00de5ab41705e509</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <child ID="175" label="ncro: bacterial_miRNA"
            created="1433778917213" x="34.0" y="68.5" width="168.0"
            height="27.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#FDE888</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/d3e2e705c0a80c4a00de5ab49ff3d652</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <child ID="176" label="ncro: virus_miRNA"
            created="1433778917214" x="34.0" y="91.75" width="140.0"
            height="27.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
            <fillColor>#FDE888</fillColor>
            <strokeColor>#000000</strokeColor>
            <textColor>#000000</textColor>
            <font>Arial-plain-16</font>
            <URIString>http://vue.tufts.edu/rdf/resource/d3e2e705c0a80c4a00de5ab442e45074</URIString>
            <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
        </child>
        <child ID="178" label="And so on..." created="1433778917214"
            x="34.0" y="115.0" width="102.0" height="18.0"
            strokeWidth="0.0" autoSized="false" xsi:type="text">
            <strokeColor>#404040</strokeColor>
            <textColor>#000000</textColor>
            <font>-plain-1</font>
            <URIString>http://vue.tufts.edu/rdf/resource/d3e2e706c0a80c4a00de5ab431fd452e</URIString>
            <richText>&lt;html&gt;
  &lt;head style="color: #000000" color="#000000"&gt;
    &lt;style type="text/css"&gt;
      &lt;!--
        body { font-family: Arial; margin-bottom: 0px; margin-right: 0px; font-size: 12; margin-top: 0px; margin-left: 0px; color: #000000 }
        ol { font-family: Arial; font-size: 12; margin-top: 6; margin-left: 30; list-style-position: outside; vertical-align: middle }
        p { margin-bottom: 0; margin-right: 0; margin-top: 0; margin-left: 0; color: #000000 }
        ul { font-family: Arial; font-size: 12; margin-top: 6; margin-left: 30; list-style-position: outside; vertical-align: middle }
      --&gt;
    &lt;/style&gt;
    
  &lt;/head&gt;
  &lt;body&gt;
    &lt;p style="color: #000000" color="#000000"&gt;
      And so on...
    &lt;/p&gt;
  &lt;/body&gt;
&lt;/html&gt;
</richText>
            <label>And so on...</label>
        </child>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="179" label="is_a" layerID="1" created="1433778987473"
        x="235.45413" y="456.25195" width="145.12277" height="56.33496"
        strokeWidth="1.0" autoSized="false" controlCount="2"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d3e52679c0a80c4a00de5ab47f2554ab</URIString>
        <point1 x="235.95413" y="497.51135"/>
        <point2 x="380.0769" y="456.75195"/>
        <ID1 xsi:type="node">172</ID1>
        <ID2 xsi:type="node">92</ID2>
        <ctrlPoint0 x="299.49146" y="508.1682" xsi:type="point"/>
        <ctrlPoint1 x="381.23645" y="521.9758" xsi:type="point"/>
    </child>
    <layer ID="1" label="Layer 1" created="1431734769931" x="0.0"
        y="0.0" width="1.4E-45" height="1.4E-45" strokeWidth="0.0" autoSized="false">
        <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b60a147a5701bb3ae4dfffb6d3</URIString>
    </layer>
    <userZoom>1.0</userZoom>
    <userOrigin x="-148.65727" y="-259.2552"/>
    <presentationBackground>#202020</presentationBackground>
    <PathwayList currentPathway="0" revealerIndex="-1">
        <pathway ID="0" label="Untitled Pathway" created="1431734769931"
            x="0.0" y="0.0" width="1.4E-45" height="1.4E-45"
            strokeWidth="0.0" autoSized="false" currentIndex="0" open="true">
            <strokeColor>#B3CC33CC</strokeColor>
            <textColor>#000000</textColor>
            <font>SansSerif-plain-14</font>
            <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b70a147a5701bb3ae4c573cfec</URIString>
            <masterSlide ID="2" created="1431734769956" x="0.0" y="0.0"
                width="800.0" height="600.0" locked="true"
                strokeWidth="0.0" autoSized="false">
                <fillColor>#000000</fillColor>
                <strokeColor>#404040</strokeColor>
                <textColor>#000000</textColor>
                <font>SansSerif-plain-14</font>
                <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b80a147a5701bb3ae4de1adb6c</URIString>
                <titleStyle ID="3" label="Header"
                    created="1431734769971" x="335.5" y="171.0"
                    width="129.0" height="58.0" strokeWidth="0.0"
                    autoSized="true" isStyle="true" xsi:type="node">
                    <strokeColor>#404040</strokeColor>
                    <textColor>#FFFFFF</textColor>
                    <font>Gill Sans-plain-36</font>
                    <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b80a147a5701bb3ae4d4c9c2ee</URIString>
                    <shape xsi:type="rectangle"/>
                </titleStyle>
                <textStyle ID="4" label="Slide Text"
                    created="1431734769972" x="346.5" y="280.5"
                    width="107.0" height="39.0" strokeWidth="0.0"
                    autoSized="true" isStyle="true" xsi:type="node">
                    <strokeColor>#404040</strokeColor>
                    <textColor>#FFFFFF</textColor>
                    <font>Gill Sans-plain-22</font>
                    <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b90a147a5701bb3ae4b6db12b2</URIString>
                    <shape xsi:type="rectangle"/>
                </textStyle>
                <linkStyle ID="5" label="Links" created="1431734769973"
                    x="373.5" y="383.0" width="53.0" height="34.0"
                    strokeWidth="0.0" autoSized="true" isStyle="true" xsi:type="node">
                    <strokeColor>#404040</strokeColor>
                    <textColor>#B3BFE3</textColor>
                    <font>Gill Sans-plain-18</font>
                    <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b90a147a5701bb3ae44e67e965</URIString>
                    <shape xsi:type="rectangle"/>
                </linkStyle>
            </masterSlide>
        </pathway>
    </PathwayList>
    <date>2015-05-15</date>
    <modelVersion>6</modelVersion>
    <saveLocation>C:\Users\Jingshan Huang\Documents\GitHub\OmniSearch\ncRO-Ontology-Design</saveLocation>
    <saveFile>C:\Users\Jingshan Huang\Documents\GitHub\OmniSearch\ncRO-Ontology-Design\ncRO-core-v5.vue</saveFile>
</LW-MAP>

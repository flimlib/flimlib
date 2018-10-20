$classifier = "natives-windows_64"

# Get version number
[xml]$pom = Get-Content 'pom.xml'
$pomNamespace = @{ pomNs = 'http://maven.apache.org/POM/4.0.0'; };
$groupId = (Select-Xml -Xml $pom -XPath '/pomNs:project/pomNs:groupId' -Namespace $pomNamespace).Node.InnerText;
$groupIdForURL = $groupId -replace "\.","/";
$artifactId = (Select-Xml -Xml $pom -XPath '/pomNs:project/pomNs:artifactId' -Namespace $pomNamespace).Node.InnerText;
$version = (Select-Xml -Xml $pom -XPath '/pomNs:project/pomNs:version' -Namespace $pomNamespace).Node.InnerText;

If (($Env:APPVEYOR_REPO_TAG -match "false")  -and ($Env:APPVEYOR_REPO_BRANCH -match "master") -and !(Test-Path Env:\APPVEYOR_PULL_REQUEST_NUMBER)) {
    "== Building and deploying master SNAPSHOT =="
    & "mvn" "-B" "-Pdeploy-to-imagej" "deploy" 2> $null
} ElseIf (($Env:APPVEYOR_REPO_TAG -match "true") -and (Test-Path ($Env:APPVEYOR_BUILD_FOLDER + "\release.properties"))) {
    "== Cutting and deploying release version =="
    & "mvn" "-B" "release:perform" 2> $null

    # Check if the parent folder in the Nexus is available
    $responseFolder = try { (Invoke-Webrequest -uri "http://maven.imagej.net/content/repositories/releases/$groupIdForURL/$artifactId/$version/" -UseBasicParsing -method head -TimeoutSec 5).statuscode } catch { $_.Exception.Response.StatusCode.Value__ }

    If ($responseFolder -ne 200) {
        exit $LASTEXITCODE
    }

    Get-ChildItem ".\target\checkout\target" -Filter *.jar |
    Foreach-Object {
        $artifactPath = $_.FullName

        $fileName = [System.IO.Path]::GetFileName("$artifactPath")

        # Skip the non-classified artifacts
        If (!("$fileName" -like "*$classifier*")) {
            return
        }
        $extension = [System.IO.Path]::GetExtension("$artifactPath") -replace "^\.", ""
        # Check if the launcher itself was already deployed
        $responseFile = try { (Invoke-Webrequest -uri "http://maven.imagej.net/content/repositories/releases/$groupIdForURL/$artifactId/$version/$fileName" -UseBasicParsing -method head -TimeoutSec 5).statuscode } catch { $_.Exception.Response.StatusCode.Value__ }
        # Deploy only iff the parent exists and the launcher does not exist
        If ($responseFile -eq 404) {
            $files="$files,$mainFile"
            $types="$types,$mainType"
            $classifiers="$classifiers,$mainClassifier"
            $mainFile="$artifactPath"
            $mainType="$extension"
            $mainClassifier="$classifier"
        }
    }
    & "mvn" "deploy:deploy-file"`
        "-Dfile=`"$mainFile`""`
        "-Dfiles=`"$files`""`
        "-DrepositoryId=`"imagej.releases`""`
        "-Durl=`"dav:https://maven.imagej.net/content/repositories/releases`""`
        "-DgeneratePom=`"false`""`
        "-DgroupId=`"$groupId`""`
        "-DartifactId=`"$artifactId`""`
        "-Dversion=`"$version`""`
        "-Dclassifier=`"$mainClassifier`""`
        "-Dclassifiers=`"$classifiers`""`
        "-Dpackaging=`"$mainType`""`
        "-Dtypes=`"$types`""
} Else {
    "== Building the artifact locally =="
    & "mvn" "-B" "-P$profiles" "install" "javadoc:aggregate-jar" 2> $null
}

exit $LASTEXITCODE

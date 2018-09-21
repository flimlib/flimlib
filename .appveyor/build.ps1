# If ($Env:PLATFORM -match "x64") {
#     $arch = "amd64"
#     $profiles = "amd64-Windows"
#     $classifier = "win64"
# } ElseIf ($Env:PLATFORM -match "x86") {
#     $arch = "x86"
#     $profiles = "!amd64-Windows,i386-Windows"
#     $classifier = "win32"
# }

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

    # # Check if the parent folder in the Nexus is available
    # $responseFolder = try { (Invoke-Webrequest -uri "http://maven.imagej.net/content/repositories/releases/$groupIdForURL/$artifactId/$version/" -UseBasicParsing -method head -TimeoutSec 5).statuscode } catch { $_.Exception.Response.StatusCode.Value__ }

    # # Check if the launcher itself was already deployed
    # $responseFile = try { (Invoke-Webrequest -uri "http://maven.imagej.net/content/repositories/releases/$groupIdForURL/$artifactId/$version/$artifactId-$version-$classifier.exe" -UseBasicParsing -method head -TimeoutSec 5).statuscode } catch { $_.Exception.Response.StatusCode.Value__ }

    # # Deploy only iff the parent exists and the launcher does not exist
    # If (($responseFolder -eq 200) -and ($responseFile -eq 404)) {
        # & "mvn" "deploy:deploy-file" "-Dfile=`"target/checkout/target/ImageJ-$classifier.exe`"" "-DrepositoryId=`"imagej.releases`"" "-Durl=`"dav:https://maven.imagej.net/content/repositories/releases`"" "-DgeneratePom=`"false`"" "-DgroupId=`"$groupId`"" "-DartifactId=`"$artifactId`"" "-Dversion=`"$version`"" "-Dclassifier=`"$classifier`"" "-Dpackaging=`"exe`"" 2> $null
    # }
} Else {
    "== Building the artifact locally =="
    & "mvn" "-B" "-P$profiles" "install" "javadoc:aggregate-jar" 2> $null
}

exit $LASTEXITCODE

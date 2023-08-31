const resource = [
    /* --- CSS --- */
    '/maxawake/assets/css/style.css',

    /* --- PWA --- */
    '/maxawake/app.js',
    '/maxawake/sw.js',

    /* --- HTML --- */
    '/maxawake/index.html',
    '/maxawake/404.html',

    
        '/maxawake/categories/',
    
        '/maxawake/tags/',
    
        '/maxawake/archives/',
    
        '/maxawake/about/',
    

    /* --- Favicons & compressed JS --- */
    
    
        '/maxawake/assets/img/favicons/android-chrome-192x192.png',
        '/maxawake/assets/img/favicons/android-chrome-512x512.png',
        '/maxawake/assets/img/favicons/apple-touch-icon.png',
        '/maxawake/assets/img/favicons/favicon-16x16.png',
        '/maxawake/assets/img/favicons/favicon-32x32.png',
        '/maxawake/assets/img/favicons/favicon.ico',
        '/maxawake/assets/img/favicons/mstile-150x150.png',
        '/maxawake/assets/js/dist/categories.min.js',
        '/maxawake/assets/js/dist/commons.min.js',
        '/maxawake/assets/js/dist/home.min.js',
        '/maxawake/assets/js/dist/misc.min.js',
        '/maxawake/assets/js/dist/page.min.js',
        '/maxawake/assets/js/dist/post.min.js'
];

/* The request url with below domain will be cached */
const allowedDomains = [
    

    'localhost:4000',

    
        'chirpy-img.netlify.app',
    

    'fonts.gstatic.com',
    'fonts.googleapis.com',
    'cdn.jsdelivr.net',
    'polyfill.io'
];

/* Requests that include the following path will be banned */
const denyUrls = [];


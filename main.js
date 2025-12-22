console.log('Javascript is loaded');


//Scroll to the right location
const scroll_motivation = document.getElementById('scroll-motivation');
const scroll_algorithm = document.getElementById('scroll-algorithm');
const scroll_usage_manual = document.getElementById('scroll-usage-manual');

scroll_motivation.addEventListener('click',()=>{
  document.querySelector('.motivation').scrollIntoView({
    behavior: "smooth"
  })
});

scroll_algorithm.addEventListener('click',()=>{
  document.querySelector('.algorithm').scrollIntoView({
    behavior: "smooth"
  })
});

scroll_usage_manual.addEventListener('click',()=>{
  document.querySelector('.usage').scrollIntoView({
    behavior: "smooth"
  })
});


//To match the size of widht of navigator items
const navigator_Items = document.querySelectorAll('.navigator-item');
// Convert NodeList to array of numbers (offsetWidth includes border)
const widths = Array.from(navigator_Items, el => el.offsetWidth);
const maxWidth = Math.max(...widths);
navigator_Items.forEach(el => {
  el.style.width = `${maxWidth}px`;
});








const heading = document.querySelector('.heading');
const closefoam = document.querySelectorAll('.closefoam');

//To make the close foam hover animation

let start = null;

;